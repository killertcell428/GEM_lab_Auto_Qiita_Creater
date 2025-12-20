"""記事スケジューラー - 週次自動投稿スケジュール管理"""
from apscheduler.schedulers.background import BackgroundScheduler
from apscheduler.triggers.cron import CronTrigger
from datetime import datetime, timedelta
from typing import Optional, Dict, Any
from pathlib import Path
import sys

# プロジェクトルートをパスに追加
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from src.config_loader import get_config
from src.crewai.state.article_state import ArticleState
from src.crewai.crews.plan_crew import execute_plan_phase
from src.crewai.crews.do_crew import execute_do_phase
from src.notifications.email_sender import EmailSender
from src.scheduler.draft_scheduler import DraftScheduler


class ArticleScheduler:
    """記事スケジューラー"""
    
    def __init__(self):
        self.scheduler = BackgroundScheduler()
        self.config = get_config()
        self.email_sender = EmailSender()
        self.draft_scheduler = DraftScheduler()
        self._setup_schedules()
    
    def _setup_schedules(self):
        """スケジュールを設定"""
        workflow_config = self.config.get("workflow", {})
        schedule_config = workflow_config.get("schedule", {})
        
        # 水曜日のスケジュール（火曜日に実行）
        wednesday_config = schedule_config.get("wednesday", {})
        if wednesday_config.get("enabled", False):
            self.scheduler.add_job(
                self._create_scheduled_article,
                CronTrigger(day_of_week="tue", hour=9, minute=0),  # 火曜日 9:00
                id="wednesday_article",
                args=[wednesday_config],
                replace_existing=True
            )
        
        # 金曜日のスケジュール（金曜日に直接投稿）
        friday_config = schedule_config.get("friday", {})
        if friday_config.get("enabled", False):
            # 金曜日の指定時刻に実行（デフォルト: 10:00）
            friday_hour = friday_config.get("hour", 10)
            friday_minute = friday_config.get("minute", 0)
            
            self.scheduler.add_job(
                self._publish_friday_draft,
                CronTrigger(day_of_week="fri", hour=friday_hour, minute=friday_minute),  # 金曜日
                id="friday_draft_article",
                replace_existing=True
            )
        
        # 承認期限チェック（毎時実行）
        self.scheduler.add_job(
            self._check_approval_deadlines,
            CronTrigger(minute=0),  # 毎時0分
            id="check_approval_deadlines",
            replace_existing=True
        )
    
    def _create_scheduled_article(self, schedule_config: Dict[str, Any]):
        """スケジュールに基づいて記事を作成"""
        try:
            article_type = schedule_config.get("article_type", "general")
            description = schedule_config.get("description", "")
            
            # トピックを生成（将来的にはAIで生成）
            topic = f"{description} - {datetime.now().strftime('%Y年%m月%d日')}"
            
            # ArticleStateを初期化
            article_id = f"article_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            article_state = ArticleState(article_id=article_id, topic=topic)
            article_state.save()
            
            # Plan Phaseを実行
            context = {
                "article_type": article_type,
                "description": description,
                "scheduled": True
            }
            article_state = execute_plan_phase(article_state, topic, context)
            
            # Do Phaseを実行
            article_state = execute_do_phase(article_state, context)
            
            # 承認待ち状態に設定
            scheduled_date = datetime.now() + timedelta(days=1)  # 翌日投稿予定
            article_state.pending_approval = True
            article_state.approval_deadline = (datetime.now() + timedelta(hours=24)).isoformat()
            article_state.approval_status = "pending"
            article_state.scheduled_publish_date = scheduled_date.isoformat()
            article_state.save()
            
            # メール通知を送信
            self._send_approval_notification(article_state)
            
            print(f"[SCHEDULER] スケジュール記事を作成しました: {article_id}")
            
        except Exception as e:
            print(f"[ERROR] スケジュール記事作成エラー: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def _send_approval_notification(self, article_state: ArticleState):
        """承認待ち通知メールを送信"""
        try:
            plan = article_state.plan if isinstance(article_state.plan, dict) else {}
            title = plan.get("title", article_state.topic or "タイトル未設定")
            
            # Web UIのURL（環境変数から取得、デフォルトはlocalhost）
            import os
            web_url = os.getenv("WEB_UI_URL", "http://localhost:3000")
            article_url = f"{web_url}/articles/{article_state.article_id}"
            
            # 承認期限
            deadline_str = "24時間以内"
            if article_state.approval_deadline:
                deadline_dt = datetime.fromisoformat(article_state.approval_deadline)
                deadline_str = deadline_dt.strftime("%Y年%m月%d日 %H:%M")
            
            # メール本文
            content_preview = article_state.content[:500] + "..." if article_state.content and len(article_state.content) > 500 else (article_state.content or "（内容準備中）")
            
            body = f"""新しい記事が作成されました。承認をお願いします。

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

【タイトル】
{title}

【承認期限】
{deadline_str}

【内容抜粋】
{content_preview}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

記事を確認・承認する: {article_url}

承認期限を過ぎると自動的に投稿されます。
"""
            
            # メール送信
            self.email_sender.send_approval_notification(
                article_title=title,
                article_url=article_url,
                approval_deadline=article_state.approval_deadline,
                article_content=article_state.content
            )
            
            print(f"[SCHEDULER] 承認通知メールを送信しました: {article_state.article_id}")
            
        except Exception as e:
            print(f"[ERROR] 承認通知メール送信エラー: {str(e)}")
    
    def _publish_friday_draft(self):
        """金曜日にドラフトファイルを直接投稿"""
        try:
            print(f"[SCHEDULER] 金曜日のドラフト投稿を開始: {datetime.now()}")
            result = self.draft_scheduler.publish_next_draft(auto_publish=True)
            
            if result.get("success"):
                print(f"[SCHEDULER] 金曜日のドラフト投稿成功: {result.get('url', 'N/A')}")
            else:
                print(f"[ERROR] 金曜日のドラフト投稿失敗: {result.get('error', 'Unknown error')}")
                
        except Exception as e:
            print(f"[ERROR] 金曜日のドラフト投稿エラー: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def _check_approval_deadlines(self):
        """承認期限をチェックして、期限超過の場合は自動投稿"""
        try:
            from src.config_loader import get_config
            from src.publisher.qiita_publisher import QiitaPublisher
            from src.generator.article_generator import generate_dynamic_tags
            
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
            
            if not articles_dir.exists():
                return
            
            now = datetime.now()
            publisher = QiitaPublisher()
            
            for json_file in articles_dir.glob("*.json"):
                try:
                    article_state = ArticleState.load(json_file.stem)
                    
                    # 承認待ちで期限超過の場合
                    if (article_state.pending_approval and 
                        article_state.approval_status == "pending" and
                        article_state.approval_deadline):
                        
                        deadline = datetime.fromisoformat(article_state.approval_deadline)
                        
                        if now > deadline:
                            # 自動投稿
                            if article_state.content:
                                plan = article_state.plan if isinstance(article_state.plan, dict) else {}
                                title = plan.get("title", article_state.topic or "タイトル未設定")
                                
                                # タグ生成
                                try:
                                    tags = generate_dynamic_tags(title, article_state.content)
                                except Exception as e:
                                    print(f"[WARN] タグ生成エラー: {str(e)}")
                                    tags = []
                                
                                article_data = {
                                    "title": title,
                                    "content": article_state.content,
                                    "tags": tags[:5] if tags else []
                                }
                                
                                result = publisher.publish_article(article_data, private=False, tweet=False)
                                
                                if result.get("success", False):
                                    article_state.qiita_url = result.get("url")
                                    article_state.qiita_item_id = result.get("id")
                                    article_state.pending_approval = False
                                    article_state.approval_status = "auto_published"
                                    article_state.save()
                                    
                                    print(f"[SCHEDULER] 承認期限超過により自動投稿しました: {article_state.article_id}")
                                else:
                                    print(f"[ERROR] 自動投稿に失敗しました: {article_state.article_id}")
                            else:
                                print(f"[WARN] 記事本文がないため自動投稿をスキップ: {article_state.article_id}")
                                
                except Exception as e:
                    print(f"[WARN] 承認期限チェックエラー ({json_file.stem}): {str(e)}")
                    continue
                    
        except Exception as e:
            print(f"[ERROR] 承認期限チェック処理エラー: {str(e)}")
            import traceback
            traceback.print_exc()
    
    def start(self):
        """スケジューラーを開始"""
        if not self.scheduler.running:
            self.scheduler.start()
            print("[SCHEDULER] スケジューラーを開始しました")
    
    def stop(self):
        """スケジューラーを停止"""
        if self.scheduler.running:
            self.scheduler.shutdown()
            print("[SCHEDULER] スケジューラーを停止しました")


# グローバルスケジューラーインスタンス
_scheduler_instance: Optional[ArticleScheduler] = None


def get_scheduler() -> ArticleScheduler:
    """スケジューラーインスタンスを取得"""
    global _scheduler_instance
    if _scheduler_instance is None:
        _scheduler_instance = ArticleScheduler()
    return _scheduler_instance

