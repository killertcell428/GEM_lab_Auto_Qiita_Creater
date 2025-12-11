"""記事サービス - CrewAI Orchestratorとの統合"""
from typing import Dict, Any, Optional, List
from pathlib import Path
import sys
import json
import traceback

# プロジェクトルートをパスに追加
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))

from src.crewai.orchestrator import PDCAOrchestrator
from src.crewai.state.article_state import ArticleState
from src.crewai.state.pdca_state import PDCAState
from src.crewai.crews.plan_crew import execute_plan_phase
from src.crewai.crews.do_crew import execute_do_phase
from src.crewai.crews.check_crew import execute_check_phase
from src.crewai.crews.act_crew import execute_act_phase
from api.app.utils.viewmodel import article_state_to_viewmodel


class ArticleService:
    """記事サービス"""
    
    def __init__(self):
        self.orchestrator = PDCAOrchestrator()
    
    def create_article(self, topic: str, draft_file: Optional[str] = None, context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """新規記事作成（Plan Phase開始）"""
        try:
            if not topic or not topic.strip():
                raise ValueError("トピックが指定されていません")
            
            if context is None:
                context = {}
            
            # ArticleStateを初期化
            from datetime import datetime
            article_id = f"article_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            article_state = ArticleState(article_id=article_id, topic=topic)
            article_state.save()
            
            # Plan Phaseを実行
            try:
                if draft_file:
                    article_state = execute_plan_phase(article_state, topic, context, draft_file)
                else:
                    article_state = execute_plan_phase(article_state, topic, context)
            except Exception as e:
                # Plan Phase実行エラーを記録
                print(f"[ERROR] Plan Phase実行エラー (article_id: {article_id}): {str(e)}")
                raise RuntimeError(f"Plan Phaseの実行に失敗しました: {str(e)}") from e
            
            if article_state is None:
                raise RuntimeError("Plan Phaseの結果が取得できませんでした（article_stateがNone）")
            # to_dict() が None を返した場合に備えた防御
            state_dict = article_state.to_dict() if hasattr(article_state, "to_dict") else None
            if state_dict is None:
                raise RuntimeError("Plan Phase結果のシリアライズに失敗しました（state_dictがNone）")
            
            # ViewModelに変換
            return article_state_to_viewmodel(state_dict)
        except ValueError as e:
            raise
        except Exception as e:
            print(f"[ERROR] 記事作成エラー: {str(e)}")
            traceback.print_exc()
            raise RuntimeError(f"記事の作成に失敗しました: {str(e)}") from e
    
    def get_article(self, article_id: str) -> Dict[str, Any]:
        """記事取得"""
        try:
            if not article_id or not article_id.strip():
                raise ValueError("article_idが指定されていません")
            
            # ファイルの存在を確認
            from src.config_loader import get_config
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
            file_path = articles_dir / f"{article_id}.json"
            
            if not file_path.exists():
                raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
            
            article_state = ArticleState.load(article_id)
            return article_state_to_viewmodel(article_state.to_dict())
        except FileNotFoundError:
            raise
        except Exception as e:
            print(f"[ERROR] 記事取得エラー (article_id: {article_id}): {str(e)}")
            raise RuntimeError(f"記事の取得に失敗しました: {str(e)}") from e
    
    def list_articles(self, limit: int = 50) -> List[Dict[str, Any]]:
        """記事一覧取得"""
        from api.app.utils.viewmodel import phase_to_status_text, get_next_action_hint
        
        articles_dir = Path("data/state/articles")
        if not articles_dir.exists():
            return []
        
        articles = []
        for json_file in sorted(articles_dir.glob("*.json"), key=lambda x: x.stat().st_mtime, reverse=True)[:limit]:
            try:
                article_state = ArticleState.load(json_file.stem)
                state_dict = article_state.to_dict()
                plan = state_dict.get("plan", {})
                
                articles.append({
                    "id": article_state.article_id,
                    "title": plan.get("title") if isinstance(plan, dict) else state_dict.get("topic", "タイトル未設定"),
                    "uiStatusText": phase_to_status_text(state_dict.get("current_phase", "plan")),
                    "phase": state_dict.get("current_phase", "plan"),
                    "nextActionHint": get_next_action_hint(state_dict.get("current_phase", "plan"), state_dict),
                    "createdAt": state_dict.get("created_at"),
                    "updatedAt": state_dict.get("updated_at"),
                    # 承認関連フィールド
                    "pendingApproval": state_dict.get("pending_approval", False),
                    "approvalDeadline": state_dict.get("approval_deadline"),
                    "approvalStatus": state_dict.get("approval_status"),
                    "scheduledPublishDate": state_dict.get("scheduled_publish_date"),
                    # KPIサマリー
                    "kpiSummary": state_dict.get("kpi")
                })
            except Exception as e:
                print(f"Warning: Failed to load article {json_file.stem}: {e}")
                continue
        
        return articles
    
    def update_article(self, article_id: str, content: Optional[str] = None, title: Optional[str] = None) -> Dict[str, Any]:
        """記事更新"""
        try:
            if not article_id or not article_id.strip():
                raise ValueError("article_idが指定されていません")
            
            # ファイルの存在を確認
            from src.config_loader import get_config
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
            file_path = articles_dir / f"{article_id}.json"
            
            if not file_path.exists():
                raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
            
            article_state = ArticleState.load(article_id)
            
            if content is not None:
                article_state.content = content
            if title is not None:
                if article_state.plan and isinstance(article_state.plan, dict):
                    article_state.plan["title"] = title
                else:
                    article_state.plan = {"title": title}
            
            article_state.save()
            return article_state_to_viewmodel(article_state.to_dict())
        except FileNotFoundError:
            raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
        except Exception as e:
            print(f"[ERROR] 記事更新エラー (article_id: {article_id}): {str(e)}")
            raise RuntimeError(f"記事の更新に失敗しました: {str(e)}") from e
    
    def delete_article(self, article_id: str) -> bool:
        """記事削除"""
        articles_dir = Path("data/state/articles")
        json_file = articles_dir / f"{article_id}.json"
        if json_file.exists():
            json_file.unlink()
            return True
        return False
    
    def execute_plan_phase(self, article_id: str, context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Plan Phase実行"""
        try:
            if not article_id or not article_id.strip():
                raise ValueError("article_idが指定されていません")
            
            # ファイルの存在を確認
            from src.config_loader import get_config
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
            file_path = articles_dir / f"{article_id}.json"
            
            if not file_path.exists():
                raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
            
            article_state = ArticleState.load(article_id)
            if context is None:
                context = {}
            
            try:
                article_state = execute_plan_phase(
                    article_state,
                    article_state.topic or "未設定",
                    context,
                    None  # draft_fileは既に読み込まれている想定
                )
            except Exception as e:
                print(f"[ERROR] Plan Phase実行エラー (article_id: {article_id}): {str(e)}")
                raise RuntimeError(f"Plan Phaseの実行に失敗しました: {str(e)}") from e
            
            return article_state_to_viewmodel(article_state.to_dict())
        except FileNotFoundError:
            raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
        except Exception as e:
            if isinstance(e, (ValueError, FileNotFoundError, RuntimeError)):
                raise
            print(f"[ERROR] Plan Phase実行エラー (article_id: {article_id}): {str(e)}")
            raise RuntimeError(f"Plan Phaseの実行に失敗しました: {str(e)}") from e
    
    def execute_do_phase(self, article_id: str, context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Do Phase実行"""
        try:
            if not article_id or not article_id.strip():
                raise ValueError("article_idが指定されていません")
            
            # ファイルの存在を確認
            from src.config_loader import get_config
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
            file_path = articles_dir / f"{article_id}.json"
            
            if not file_path.exists():
                raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
            
            article_state = ArticleState.load(article_id)
            if context is None:
                context = {}
            
            try:
                article_state = execute_do_phase(article_state, context)
            except Exception as e:
                print(f"[ERROR] Do Phase実行エラー (article_id: {article_id}): {str(e)}")
                raise RuntimeError(f"Do Phaseの実行に失敗しました: {str(e)}") from e
            
            return article_state_to_viewmodel(article_state.to_dict())
        except FileNotFoundError:
            raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
        except Exception as e:
            if isinstance(e, (ValueError, FileNotFoundError, RuntimeError)):
                raise
            print(f"[ERROR] Do Phase実行エラー (article_id: {article_id}): {str(e)}")
            raise RuntimeError(f"Do Phaseの実行に失敗しました: {str(e)}") from e
    
    def execute_check_phase(self, article_id: str, context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Check Phase実行"""
        try:
            if not article_id or not article_id.strip():
                raise ValueError("article_idが指定されていません")
            
            # ファイルの存在を確認
            from src.config_loader import get_config
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
            file_path = articles_dir / f"{article_id}.json"
            
            if not file_path.exists():
                raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
            
            article_state = ArticleState.load(article_id)
            if context is None:
                context = {}
            
            try:
                article_state = execute_check_phase(article_state, context)
            except Exception as e:
                print(f"[ERROR] Check Phase実行エラー (article_id: {article_id}): {str(e)}")
                raise RuntimeError(f"Check Phaseの実行に失敗しました: {str(e)}") from e
            
            return article_state_to_viewmodel(article_state.to_dict())
        except FileNotFoundError:
            raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
        except Exception as e:
            if isinstance(e, (ValueError, FileNotFoundError, RuntimeError)):
                raise
            print(f"[ERROR] Check Phase実行エラー (article_id: {article_id}): {str(e)}")
            raise RuntimeError(f"Check Phaseの実行に失敗しました: {str(e)}") from e
    
    def execute_act_phase(self, article_id: str, context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Act Phase実行"""
        try:
            if not article_id or not article_id.strip():
                raise ValueError("article_idが指定されていません")
            
            # ファイルの存在を確認
            from src.config_loader import get_config
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
            file_path = articles_dir / f"{article_id}.json"
            
            if not file_path.exists():
                raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
            
            article_state = ArticleState.load(article_id)
            
            # article_idからPDCAStateを検索
            from src.config_loader import get_config
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            pdca_dir = Path(state_config.get("pdca_dir", "data/state/pdca"))
            
            pdca_state = None
            if pdca_dir.exists():
                for json_file in pdca_dir.glob("*.json"):
                    try:
                        with open(json_file, "r", encoding="utf-8") as f:
                            data = json.load(f)
                            if data.get("article_id") == article_id:
                                pdca_state = PDCAState(**data)
                                break
                    except Exception as e:
                        print(f"[WARN] PDCAState読み込みエラー ({json_file}): {str(e)}")
                        continue
            
            # PDCAStateが見つからない場合は新規作成
            if pdca_state is None:
                from datetime import datetime
                cycle_id = f"cycle_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                pdca_state = PDCAState(cycle_id=cycle_id, article_id=article_id)
            
            if context is None:
                context = {}
            
            try:
                article_state, pdca_state = execute_act_phase(article_state, pdca_state, context)
            except Exception as e:
                print(f"[ERROR] Act Phase実行エラー (article_id: {article_id}): {str(e)}")
                raise RuntimeError(f"Act Phaseの実行に失敗しました: {str(e)}") from e
            
            return article_state_to_viewmodel(article_state.to_dict())
        except FileNotFoundError:
            raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
        except Exception as e:
            if isinstance(e, (ValueError, FileNotFoundError, RuntimeError)):
                raise
            print(f"[ERROR] Act Phase実行エラー (article_id: {article_id}): {str(e)}")
            raise RuntimeError(f"Act Phaseの実行に失敗しました: {str(e)}") from e
    
    def execute_full_pdca(self, article_id: str, auto_publish: bool = False, context: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """フルPDCAサイクル実行"""
        try:
            if not article_id or not article_id.strip():
                raise ValueError("article_idが指定されていません")
            
            # ファイルの存在を確認
            from src.config_loader import get_config
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
            file_path = articles_dir / f"{article_id}.json"
            
            if not file_path.exists():
                raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
            
            article_state = ArticleState.load(article_id)
            if context is None:
                context = {}
            
            try:
                result = self.orchestrator.execute_full_pdca_cycle(
                    article_state.topic or "未設定",
                    None,  # draft_file
                    auto_publish,
                    context
                )
            except Exception as e:
                print(f"[ERROR] フルPDCAサイクル実行エラー (article_id: {article_id}): {str(e)}")
                raise RuntimeError(f"フルPDCAサイクルの実行に失敗しました: {str(e)}") from e
            
            # 最新の状態を取得
            article_state = ArticleState.load(article_id)
            return article_state_to_viewmodel(article_state.to_dict())
        except FileNotFoundError:
            raise
        except Exception as e:
            if isinstance(e, (ValueError, FileNotFoundError, RuntimeError)):
                raise
            print(f"[ERROR] フルPDCAサイクル実行エラー (article_id: {article_id}): {str(e)}")
            raise RuntimeError(f"フルPDCAサイクルの実行に失敗しました: {str(e)}") from e
    
    def publish_to_qiita(self, article_id: str, private: bool = False, tweet: bool = False) -> Dict[str, Any]:
        """Qiita投稿"""
        try:
            if not article_id or not article_id.strip():
                raise ValueError("article_idが指定されていません")
            
            # ファイルの存在を確認
            from src.config_loader import get_config
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
            file_path = articles_dir / f"{article_id}.json"
            
            if not file_path.exists():
                raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
            
            article_state = ArticleState.load(article_id)
            
            if not article_state.content:
                raise ValueError("記事本文がありません")
            
            from src.publisher.qiita_publisher import QiitaPublisher
            from src.generator.article_generator import generate_dynamic_tags
            
            publisher = QiitaPublisher()
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
            
            result = publisher.publish_article(article_data, private=private, tweet=tweet)
            
            if not result.get("success", False):
                error_msg = result.get("error", "不明なエラー")
                raise RuntimeError(f"Qiitaへの投稿に失敗しました: {error_msg}")
            
            # ArticleStateを更新
            article_state.qiita_url = result.get("url")
            article_state.qiita_item_id = result.get("id")
            article_state.save()
            
            return {
                "success": True,
                "url": result.get("url"),
                "id": result.get("id"),
                "message": "Qiitaに投稿しました"
            }
        except FileNotFoundError:
            raise FileNotFoundError(f"記事が見つかりません (article_id: {article_id})")
        except (ValueError, RuntimeError) as e:
            raise
        except Exception as e:
            print(f"[ERROR] Qiita投稿エラー (article_id: {article_id}): {str(e)}")
            raise RuntimeError(f"Qiitaへの投稿に失敗しました: {str(e)}") from e

