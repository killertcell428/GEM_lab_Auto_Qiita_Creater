"""記事関連のPydanticモデル"""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime


class ArticleCreateRequest(BaseModel):
    """新規記事作成リクエスト"""
    topic: str = Field(..., description="記事のトピック")
    draft_file: Optional[str] = Field(None, description="ドラフトファイルのパス")
    target_audience: Optional[str] = Field(None, description="対象読者")
    article_tone: Optional[str] = Field(None, description="記事トーン")
    code_ratio: Optional[str] = Field(None, description="コード比率")
    theory_depth: Optional[str] = Field(None, description="理論の深さ")
    environment: Optional[str] = Field(None, description="想定環境")


class ArticleUpdateRequest(BaseModel):
    """記事更新リクエスト"""
    content: Optional[str] = Field(None, description="記事本文（Markdown）")
    title: Optional[str] = Field(None, description="記事タイトル")

class FeedbackCreateRequest(BaseModel):
    """フィードバック作成リクエスト"""
    content: str
    target_section: Optional[str] = Field(None, description="対象セクション")
    intent: Optional[str] = Field(None, description="意図")
    priority: Optional[int] = Field(None, description="優先度（0-100）")
    phase: Optional[str] = Field(None, description="対象フェーズ (plan/do/check/act/publish/analyze/other)")

class FeedbackResponse(BaseModel):
    """フィードバックレスポンス"""
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


class ArticleViewModel(BaseModel):
    """記事ViewModel（UI用）"""
    id: str
    title: str
    uiStatusText: str
    phase: str
    markdown: str
    plan: Optional[Dict[str, Any]] = None
    reviewResult: Optional[Dict[str, Any]] = None
    qiitaUrl: Optional[str] = None
    qiitaItemId: Optional[str] = None
    kpiSummary: Optional[Dict[str, Any]] = None
    analysisResults: Optional[Dict[str, Any]] = None
    nextActionHint: str
    createdAt: Optional[str] = None
    updatedAt: Optional[str] = None
    feedbackHistory: List[Dict[str, Any]] = []
    # 承認関連フィールド
    pendingApproval: bool = False
    approvalDeadline: Optional[str] = None
    approvalStatus: Optional[str] = None
    scheduledPublishDate: Optional[str] = None
    # リサーチ結果
    researchReport: Optional[str] = None


class ArticleListItem(BaseModel):
    """記事一覧アイテム"""
    id: str
    title: str
    uiStatusText: str
    phase: str
    nextActionHint: str
    createdAt: Optional[str] = None
    updatedAt: Optional[str] = None

