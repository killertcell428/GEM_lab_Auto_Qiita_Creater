"""Human Feedback APIルーター"""
from fastapi import APIRouter, HTTPException
from typing import List
from api.app.models.feedback import FeedbackCreateRequest, FeedbackResponse
from src.crewai.human_loop import HumanLoop
from src.crewai.state.article_state import ArticleState
from src.crewai.state.feedback_queue import HumanFeedback
from datetime import datetime

router = APIRouter()
human_loop = HumanLoop()


@router.post("/{article_id}/feedback", response_model=FeedbackResponse, status_code=201)
async def create_feedback(article_id: str, request: FeedbackCreateRequest):
    """Human Feedback追加"""
    try:
        # バリデーション
        if not request.content or not request.content.strip():
            raise HTTPException(status_code=400, detail="フィードバック内容を入力してください")
        
        article_state = ArticleState.load(article_id)
        phase = article_state.current_phase
        
        feedback = human_loop.add_feedback(
            article_id=article_id,
            phase=phase,
            content=request.content.strip(),
            feedback_type="instruction",
            target_agent=None,  # Human Feedback Agentが判定
            priority=request.priority
        )
        # ArticleStateにも即時保存して履歴に反映
        article_state.add_human_feedback({
            "feedback_id": feedback.feedback_id,
            "content": request.content.strip(),
            "target_section": request.target_section,
            "intent": request.intent,
            "priority": request.priority or 5,
            "status": "pending",
            "phase": phase,
            "created_at": datetime.now().isoformat(),
            "processed_at": None
        })
        article_state.save()
        
        # フィードバックを処理（バックグラウンドで実行される想定）
        # 実際の処理は非同期で実行される
        
        return FeedbackResponse(
            feedback_id=feedback.feedback_id,
            article_id=feedback.article_id,
            content=feedback.content,
            target_section=request.target_section,
            intent=request.intent,
            priority=feedback.priority,
            status=feedback.status,
            created_at=feedback.created_at,
            processed_at=feedback.processed_at
        )
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"フィードバックの追加に失敗しました: {str(e)}")


@router.get("/{article_id}/feedback", response_model=List[FeedbackResponse])
async def get_feedback_history(article_id: str):
    """Feedback履歴取得"""
    try:
        article_state = ArticleState.load(article_id)
        feedback_history = article_state.human_feedback or []
        
        return [
            FeedbackResponse(
                feedback_id=fb.get("feedback_id", ""),
                article_id=article_id,
                content=fb.get("content", ""),
                target_section=fb.get("target_section"),
                intent=fb.get("intent"),
                priority=fb.get("priority", 5),
                status=fb.get("status", "pending"),
                created_at=fb.get("created_at", ""),
                processed_at=fb.get("processed_at")
            )
            for fb in feedback_history
        ]
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.put("/{article_id}/feedback/{feedback_id}", response_model=FeedbackResponse)
async def update_feedback(article_id: str, feedback_id: str, request: FeedbackCreateRequest):
    """Feedback更新"""
    try:
        article_state = ArticleState.load(article_id)
        feedback_history = article_state.human_feedback or []
        
        # フィードバックを検索して更新
        updated = False
        for fb in feedback_history:
            if fb.get("feedback_id") == feedback_id:
                fb["content"] = request.content
                if request.target_section:
                    fb["target_section"] = request.target_section
                if request.intent:
                    fb["intent"] = request.intent
                if request.priority:
                    fb["priority"] = request.priority
                updated = True
                break
        
        if not updated:
            raise HTTPException(status_code=404, detail="フィードバックが見つかりません")
        
        article_state.human_feedback = feedback_history
        article_state.save()
        
        # 更新されたフィードバックを返す
        updated_fb = next((fb for fb in feedback_history if fb.get("feedback_id") == feedback_id), None)
        return FeedbackResponse(
            feedback_id=updated_fb.get("feedback_id", ""),
            article_id=article_id,
            content=updated_fb.get("content", ""),
            target_section=updated_fb.get("target_section"),
            intent=updated_fb.get("intent"),
            priority=updated_fb.get("priority", 5),
            status=updated_fb.get("status", "pending"),
            created_at=updated_fb.get("created_at", ""),
            processed_at=updated_fb.get("processed_at")
        )
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail="記事が見つかりません")
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
