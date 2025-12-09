"""ストリーミングAPIルーター"""
from fastapi import APIRouter
from sse_starlette.sse import EventSourceResponse
from typing import AsyncGenerator
import asyncio
import json

router = APIRouter()


@router.get("/{article_id}/stream")
async def stream_article_events(article_id: str):
    """SSEでAgent実行状況をストリーミング"""
    async def event_generator() -> AsyncGenerator[str, None]:
        """イベント生成器"""
        # 実際の実装では、Agent実行中のログを監視してストリーミング
        # ここでは簡易的な実装
        
        yield f"data: {json.dumps({'type': 'connected', 'article_id': article_id})}\n\n"
        
        # 定期的に状態をチェックして送信
        from src.crewai.state.article_state import ArticleState
        from api.app.utils.viewmodel import phase_to_status_text
        
        last_phase = None
        while True:
            try:
                article_state = ArticleState.load(article_id)
                current_phase = article_state.current_phase
                
                if current_phase != last_phase:
                    yield f"data: {json.dumps({
                        'type': 'phase_change',
                        'phase': current_phase,
                        'status_text': phase_to_status_text(current_phase)
                    })}\n\n"
                    last_phase = current_phase
                
                await asyncio.sleep(2)  # 2秒ごとにチェック
            except FileNotFoundError:
                yield f"data: {json.dumps({'type': 'error', 'message': '記事が見つかりません'})}\n\n"
                break
            except Exception as e:
                yield f"data: {json.dumps({'type': 'error', 'message': str(e)})}\n\n"
                break
    
    return EventSourceResponse(event_generator())

