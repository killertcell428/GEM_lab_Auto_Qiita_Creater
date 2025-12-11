"""承認フローAPIルーター"""
from fastapi import APIRouter, HTTPException
from typing import List
from api.app.models.article import ArticleListItem
from api.app.services.article_service import ArticleService
from datetime import datetime, timedelta
from src.crewai.state.article_state import ArticleState

router = APIRouter()
service = ArticleService()


@router.get("/pending-approval", response_model=List[ArticleListItem])
async def get_pending_approval_articles():
    """承認待ち記事一覧取得"""
    try:
        from api.app.utils.viewmodel import phase_to_status_text, get_next_action_hint
        from pathlib import Path
        
        articles_dir = Path("data/state/articles")
        if not articles_dir.exists():
            return []
        
        pending_articles = []
        for json_file in sorted(articles_dir.glob("*.json"), key=lambda x: x.stat().st_mtime, reverse=True):
            try:
                article_state = ArticleState.load(json_file.stem)
                if article_state.pending_approval and article_state.approval_status == "pending":
                    state_dict = article_state.to_dict()
                    plan = state_dict.get("plan", {})
                    
                    pending_articles.append({
                        "id": article_state.article_id,
                        "title": plan.get("title") if isinstance(plan, dict) else state_dict.get("topic", "タイトル未設定"),
                        "uiStatusText": phase_to_status_text(state_dict.get("current_phase", "do")),
                        "phase": state_dict.get("current_phase", "do"),
                        "nextActionHint": f"承認期限: {article_state.approval_deadline}" if article_state.approval_deadline else "承認待ち",
                        "createdAt": state_dict.get("created_at"),
                        "updatedAt": state_dict.get("updated_at"),
                        # 承認関連フィールド
                        "pendingApproval": True,
                        "approvalDeadline": article_state.approval_deadline,
                        "approvalStatus": article_state.approval_status,
                        "scheduledPublishDate": article_state.scheduled_publish_date
                    })
            except Exception as e:
                print(f"Warning: Failed to load article {json_file.stem}: {e}")
                continue
        
        return pending_articles
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{article_id}/approve")
async def approve_article(article_id: str):
    """記事承認"""
    try:
        from src.config_loader import get_config
        from pathlib import Path
        from src.publisher.qiita_publisher import QiitaPublisher
        from src.generator.article_generator import generate_dynamic_tags
        
        config = get_config()
        crewai_config = config.get("crewai", {})
        state_config = crewai_config.get("state", {})
        articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
        file_path = articles_dir / f"{article_id}.json"
        
        if not file_path.exists():
            raise HTTPException(status_code=404, detail="記事が見つかりません")
        
        article_state = ArticleState.load(article_id)
        
        if not article_state.pending_approval:
            raise HTTPException(status_code=400, detail="この記事は承認待ちではありません")
        
        if article_state.approval_status != "pending":
            raise HTTPException(status_code=400, detail=f"この記事は既に{article_state.approval_status}です")
        
        # 承認処理: Qiitaに投稿
        if not article_state.content:
            raise HTTPException(status_code=400, detail="記事本文がありません")
        
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
        
        result = publisher.publish_article(article_data, private=False, tweet=False)
        
        if not result.get("success", False):
            error_msg = result.get("error", "不明なエラー")
            raise HTTPException(status_code=500, detail=f"Qiitaへの投稿に失敗しました: {error_msg}")
        
        # ArticleStateを更新
        article_state.qiita_url = result.get("url")
        article_state.qiita_item_id = result.get("id")
        article_state.pending_approval = False
        article_state.approval_status = "approved"
        article_state.save()
        
        return {
            "success": True,
            "url": result.get("url"),
            "id": result.get("id"),
            "message": "記事を承認し、Qiitaに投稿しました"
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"承認処理に失敗しました: {str(e)}")

