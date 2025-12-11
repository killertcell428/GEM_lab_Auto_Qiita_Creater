"""Human Feedback関連のPydanticモデル"""
from pydantic import BaseModel, Field
from typing import Optional, Literal


class FeedbackCreateRequest(BaseModel):
    """Human Feedback作成リクエスト"""
    content: str = Field(..., description="フィードバック内容（自然文）")
    target_section: Optional[str] = Field(None, description="対象セクション（全体/このセクション/AIに任せる）")
    intent: Optional[Literal["修正したい", "もっと詳しく", "方針を変えたい"]] = Field(None, description="意図")
    priority: int = Field(5, ge=1, le=10, description="優先度（1-10、10が最高）")
    phase: Optional[str] = Field(None, description="対象フェーズ(plan/do/check/act/publish/analyze/other)")


class FeedbackResponse(BaseModel):
    """Human Feedbackレスポンス"""
    feedback_id: str
    article_id: str
    content: str
    target_section: Optional[str] = None
    intent: Optional[str] = None
    priority: int
    status: str
    created_at: str
    processed_at: Optional[str] = None
    phase: Optional[str] = None

