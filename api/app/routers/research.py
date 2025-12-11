"""リサーチ結果APIルーター"""
from fastapi import APIRouter, HTTPException
from typing import List, Dict, Any
from pathlib import Path
from src.crewai.state.article_state import ArticleState

router = APIRouter()


@router.get("/{article_id}/research")
async def get_article_research(article_id: str) -> Dict[str, Any]:
    """個別記事のリサーチ結果取得"""
    try:
        from src.config_loader import get_config
        
        config = get_config()
        crewai_config = config.get("crewai", {})
        state_config = crewai_config.get("state", {})
        articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
        file_path = articles_dir / f"{article_id}.json"
        
        if not file_path.exists():
            raise HTTPException(status_code=404, detail="記事が見つかりません")
        
        article_state = ArticleState.load(article_id)
        
        return {
            "article_id": article_id,
            "topic": article_state.topic,
            "research_report": article_state.research_report,
            "plan": article_state.plan,
            "created_at": article_state.created_at,
            "updated_at": article_state.updated_at
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"リサーチ結果の取得に失敗しました: {str(e)}")


@router.get("/research-library")
async def get_research_library() -> List[Dict[str, Any]]:
    """全記事のリサーチ結果一覧取得"""
    try:
        from src.config_loader import get_config
        
        config = get_config()
        crewai_config = config.get("crewai", {})
        state_config = crewai_config.get("state", {})
        articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
        
        if not articles_dir.exists():
            return []
        
        research_list = []
        for json_file in sorted(articles_dir.glob("*.json"), key=lambda x: x.stat().st_mtime, reverse=True):
            try:
                article_state = ArticleState.load(json_file.stem)
                
                # リサーチ結果がある記事のみ追加
                if article_state.research_report or article_state.plan:
                    plan = article_state.plan if isinstance(article_state.plan, dict) else {}
                    
                    research_list.append({
                        "article_id": article_state.article_id,
                        "title": plan.get("title") if isinstance(plan, dict) else article_state.topic or "タイトル未設定",
                        "topic": article_state.topic,
                        "research_report": article_state.research_report,
                        "has_research": bool(article_state.research_report),
                        "has_plan": bool(article_state.plan),
                        "created_at": article_state.created_at,
                        "updated_at": article_state.updated_at
                    })
            except Exception as e:
                print(f"Warning: Failed to load article {json_file.stem}: {e}")
                continue
        
        return research_list
    except Exception as e:
        print(f"[ERROR] リサーチライブラリの取得に失敗しました: {str(e)}")
        return []  # エラー時も空配列を返す

