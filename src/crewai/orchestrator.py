"""メインオーケストレーター - 全Phaseを統合したPDCAサイクル"""
from typing import Dict, Any, Optional, List
from datetime import datetime
from src.crewai.crews.plan_crew import execute_plan_phase
from src.crewai.crews.do_crew import execute_do_phase
from src.crewai.crews.check_crew import execute_check_phase
from src.crewai.crews.act_crew import execute_act_phase
from src.crewai.state.article_state import ArticleState
from src.crewai.state.pdca_state import PDCAState
from src.crewai.state.feedback_queue import HumanFeedback
from src.crewai.human_loop import HumanLoop
from src.publisher.qiita_publisher import QiitaPublisher
from src.generator.article_generator import generate_dynamic_tags
from src.config_loader import get_config
import uuid


class PDCAOrchestrator:
    """PDCAサイクルを管理するオーケストレーター"""
    
    def __init__(self):
        """初期化"""
        self.human_loop = HumanLoop()
        self.publisher = QiitaPublisher()
    
    def execute_full_pdca_cycle(
        self,
        topic: str,
        draft_file: Optional[str] = None,
        auto_publish: bool = False,
        context: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        フルPDCAサイクルを実行
        
        Args:
            topic: 記事のトピック
            draft_file: ドラフトファイル（オプション、Plan Phaseで使用）
            auto_publish: 自動投稿するかどうか
            context: 追加のコンテキスト情報
        
        Returns:
            Dict[str, Any]: 実行結果
        """
        if context is None:
            context = {}
        
        # Article StateとPDCA Stateを初期化
        article_id = f"article_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        article_state = ArticleState(article_id=article_id)
        cycle_id = f"cycle_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        pdca_state = PDCAState(cycle_id=cycle_id, article_id=article_id)
        
        result = {
            "success": False,
            "article_id": article_id,
            "cycle_id": cycle_id,
            "phases": {},
            "error": None
        }
        
        try:
            # Plan Phase
            print("\n" + "="*50)
            print("[PDCA] Plan Phase開始")
            print("="*50)
            article_state = execute_plan_phase(article_state, topic, context, draft_file=draft_file)
            pdca_state.update_plan({
                "topic": topic,
                "research_report": article_state.research_report,
                "plan": article_state.plan
            })
            pdca_state.save()
            result["phases"]["plan"] = {"success": True}
            
            # Human Feedbackチェック
            pending_feedbacks = self.human_loop.check_and_process_pending_feedbacks(article_state)
            if pending_feedbacks:
                print(f"[HUMAN] {len(pending_feedbacks)}件のフィードバックを処理しました")
            
            # Do Phase
            print("\n" + "="*50)
            print("[PDCA] Do Phase開始")
            print("="*50)
            article_state = execute_do_phase(article_state, context)
            pdca_state.update_do({
                "content": article_state.content,
                "review_result": article_state.review_result
            })
            pdca_state.save()
            result["phases"]["do"] = {"success": True}
            
            # Human Feedbackチェック
            pending_feedbacks = self.human_loop.check_and_process_pending_feedbacks(article_state)
            if pending_feedbacks:
                print(f"[HUMAN] {len(pending_feedbacks)}件のフィードバックを処理しました")
            
            # 投稿（オプション）
            if auto_publish:
                print("\n" + "="*50)
                print("[PDCA] Qiitaへの投稿")
                print("="*50)
                publish_result = self._publish_article(article_state)
                if publish_result.get("success"):
                    article_state.qiita_url = publish_result.get("url")
                    article_state.qiita_item_id = publish_result.get("item_id")
                    article_state.save()
                    result["phases"]["publish"] = {"success": True, "url": publish_result.get("url")}
                else:
                    result["phases"]["publish"] = {"success": False, "error": publish_result.get("error")}
            
            # Check Phase（投稿済みの場合のみ）
            if article_state.qiita_url:
                print("\n" + "="*50)
                print("[PDCA] Check Phase開始")
                print("="*50)
                article_state = execute_check_phase(article_state, context)
                pdca_state.update_check({
                    "analysis_results": article_state.analysis_results
                })
                pdca_state.save()
                result["phases"]["check"] = {"success": True}
            else:
                print("[PDCA] Check Phaseをスキップ（記事が未投稿）")
                result["phases"]["check"] = {"skipped": True}
            
            # Act Phase
            print("\n" + "="*50)
            print("[PDCA] Act Phase開始")
            print("="*50)
            article_state, pdca_state = execute_act_phase(article_state, pdca_state, context)
            result["phases"]["act"] = {"success": True}
            
            result["success"] = True
            result["article_state"] = article_state.to_dict()
            result["pdca_state"] = pdca_state.to_dict()
            
            print("\n" + "="*50)
            print("[PDCA] 全Phase完了")
            print("="*50)
            
        except Exception as e:
            result["error"] = str(e)
            result["success"] = False
            print(f"[ERROR] PDCAサイクル実行中にエラーが発生しました: {str(e)}")
            import traceback
            traceback.print_exc()
        
        return result
    
    def _publish_article(self, article_state: ArticleState) -> Dict[str, Any]:
        """
        記事をQiitaに投稿
        
        Args:
            article_state: Article State
        
        Returns:
            Dict[str, Any]: 投稿結果
        """
        try:
            # タイトルとタグを取得
            plan = article_state.plan or {}
            title = plan.get("title", "タイトル未設定")
            if article_state.content and article_state.content.startswith("# "):
                title = article_state.content.split("\n")[0].replace("# ", "").strip()
            
            # タグを生成
            config = get_config()
            article_config = config.get("article", {})
            
            # デフォルトタグ
            default_tags = article_config.get("default_tags", ["Python", "自動化"])
            
            # planから取得したタグ
            plan_tags = plan.get("tags", [])
            if isinstance(plan_tags, str):
                plan_tags = [tag.strip() for tag in plan_tags.split(",") if tag.strip()]
            
            # 動的タグ生成（設定で有効な場合）
            all_tags = list(set(default_tags + plan_tags))
            
            if article_config.get("dynamic_tags", True) and article_state.content:
                try:
                    print("[TAG] 動的タグを生成中...")
                    dynamic_tags = generate_dynamic_tags(title, article_state.content)
                    if dynamic_tags:
                        print(f"[TAG] 生成されたタグ: {dynamic_tags}")
                        all_tags = list(set(all_tags + dynamic_tags))
                    else:
                        print("[TAG] 動的タグの生成に失敗しました（空のリスト）")
                except Exception as e:
                    print(f"[WARN] 動的タグの生成に失敗しました: {str(e)}")
            
            # 最大タグ数までに制限（Qiita APIの制限は5つ）
            max_tags = min(article_config.get("max_tags", 10), 5)  # Qiita APIの制限
            tags = all_tags[:max_tags]
            
            print(f"[TAG] 最終タグ: {tags}")
            
            # 投稿
            article_data = {
                "title": title,
                "content": article_state.content or "",
                "tags": tags
            }
            
            publish_result = self.publisher.publish_article(article_data)
            return publish_result
            
        except Exception as e:
            return {
                "success": False,
                "error": str(e)
            }
    
    def add_human_feedback(
        self,
        article_id: str,
        phase: str,
        content: str,
        feedback_type: str = "instruction",
        priority: int = 5
    ) -> HumanFeedback:
        """
        人間のフィードバックを追加
        
        Args:
            article_id: 記事ID
            phase: 現在のPhase
            content: フィードバック内容
            feedback_type: フィードバックタイプ
            priority: 優先度
        
        Returns:
            HumanFeedback: 追加されたフィードバック
        """
        return self.human_loop.add_feedback(
            article_id=article_id,
            phase=phase,
            content=content,
            feedback_type=feedback_type,
            priority=priority
        )

