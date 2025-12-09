"""Phase実行APIルーター"""
from fastapi import APIRouter, HTTPException, BackgroundTasks
from typing import Dict, Any, Optional
from api.app.models.article import ArticleViewModel
from api.app.services.article_service import ArticleService
import asyncio
from concurrent.futures import ThreadPoolExecutor

router = APIRouter()
service = ArticleService()
executor = ThreadPoolExecutor(max_workers=4)


def run_phase_sync(phase_func, *args, **kwargs):
    """同期関数を実行するヘルパー"""
    return phase_func(*args, **kwargs)


@router.post("/{article_id}/plan", response_model=ArticleViewModel)
async def execute_plan(article_id: str, background_tasks: BackgroundTasks, context: Optional[Dict[str, Any]] = None):
    """Plan Phase実行"""
    try:
        # 非同期で実行（長時間処理のため）
        loop = asyncio.get_event_loop()
        article = await loop.run_in_executor(
            executor,
            run_phase_sync,
            service.execute_plan_phase,
            article_id,
            context
        )
        return article
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{article_id}/do", response_model=ArticleViewModel)
async def execute_do(article_id: str, background_tasks: BackgroundTasks, context: Optional[Dict[str, Any]] = None):
    """Do Phase実行"""
    try:
        loop = asyncio.get_event_loop()
        article = await loop.run_in_executor(
            executor,
            run_phase_sync,
            service.execute_do_phase,
            article_id,
            context
        )
        return article
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{article_id}/check", response_model=ArticleViewModel)
async def execute_check(article_id: str, background_tasks: BackgroundTasks, context: Optional[Dict[str, Any]] = None):
    """Check Phase実行"""
    try:
        loop = asyncio.get_event_loop()
        article = await loop.run_in_executor(
            executor,
            run_phase_sync,
            service.execute_check_phase,
            article_id,
            context
        )
        return article
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{article_id}/act", response_model=ArticleViewModel)
async def execute_act(article_id: str, background_tasks: BackgroundTasks, context: Optional[Dict[str, Any]] = None):
    """Act Phase実行"""
    try:
        loop = asyncio.get_event_loop()
        article = await loop.run_in_executor(
            executor,
            run_phase_sync,
            service.execute_act_phase,
            article_id,
            context
        )
        return article
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{article_id}/pdca", response_model=ArticleViewModel)
async def execute_pdca(article_id: str, background_tasks: BackgroundTasks, auto_publish: bool = False, context: Optional[Dict[str, Any]] = None):
    """フルPDCAサイクル実行"""
    try:
        loop = asyncio.get_event_loop()
        article = await loop.run_in_executor(
            executor,
            run_phase_sync,
            service.execute_full_pdca,
            article_id,
            auto_publish,
            context
        )
        return article
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{article_id}/analyze")
async def analyze_article(article_id: str, background_tasks: BackgroundTasks):
    """分析実行（Check Phaseのエイリアス）"""
    try:
        loop = asyncio.get_event_loop()
        article = await loop.run_in_executor(
            executor,
            run_phase_sync,
            service.execute_check_phase,
            article_id,
            None
        )
        return article
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{article_id}/kpi")
async def get_kpi(article_id: str):
    """KPI取得"""
    try:
        article = service.get_article(article_id)
        return {
            "kpi": article.get("kpiSummary"),
            "analysisResults": article.get("analysisResults")
        }
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

