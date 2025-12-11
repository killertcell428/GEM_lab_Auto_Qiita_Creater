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

