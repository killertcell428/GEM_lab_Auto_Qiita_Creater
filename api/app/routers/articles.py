"""記事管理APIルーター"""
from fastapi import APIRouter, HTTPException, BackgroundTasks
from typing import List, Optional
from api.app.models.article import ArticleCreateRequest, ArticleUpdateRequest, ArticleViewModel, ArticleListItem
from api.app.services.article_service import ArticleService
import asyncio
from concurrent.futures import ThreadPoolExecutor

router = APIRouter()
service = ArticleService()
executor = ThreadPoolExecutor(max_workers=4)


def run_service_sync(service_func, *args, **kwargs):
    """同期関数を実行するヘルパー"""
    return service_func(*args, **kwargs)


@router.post("", response_model=ArticleViewModel, status_code=201)
async def create_article(request: ArticleCreateRequest, background_tasks: BackgroundTasks):
    """新規記事作成（Plan Phase開始）"""
    try:
        # バリデーション
        if not request.topic or not request.topic.strip():
            raise HTTPException(status_code=400, detail="記事のテーマを入力してください")
        
        context = {
            "target_audience": request.target_audience,
            "article_tone": request.article_tone,
            "code_ratio": request.code_ratio,
            "theory_depth": request.theory_depth,
            "environment": request.environment
        }
        
        # 非同期でPlan Phaseを実行（長時間処理のため）
        loop = asyncio.get_event_loop()
        article = await loop.run_in_executor(
            executor,
            run_service_sync,
            service.create_article,
            request.topic.strip(),
            request.draft_file,
            context
        )
        return article
    except HTTPException:
        raise
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"記事の作成に失敗しました: {str(e)}")


@router.get("", response_model=List[ArticleListItem])
async def list_articles(limit: int = 50):
    """記事一覧取得"""
    try:
        articles = service.list_articles(limit=limit)
        return articles
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{article_id}", response_model=ArticleViewModel)
async def get_article(article_id: str):
    """記事詳細取得"""
    try:
        article = service.get_article(article_id)
        return article
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.put("/{article_id}", response_model=ArticleViewModel)
async def update_article(article_id: str, request: ArticleUpdateRequest):
    """記事更新"""
    try:
        article = service.update_article(
            article_id=article_id,
            content=request.content,
            title=request.title
        )
        return article
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.delete("/{article_id}", status_code=204)
async def delete_article(article_id: str):
    """記事削除"""
    try:
        success = service.delete_article(article_id)
        if not success:
            raise HTTPException(status_code=404, detail="記事が見つかりません")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{article_id}/publish")
async def publish_article(article_id: str, private: bool = False, tweet: bool = False):
    """Qiita投稿"""
    try:
        result = service.publish_to_qiita(article_id, private=private, tweet=tweet)
        return result
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

