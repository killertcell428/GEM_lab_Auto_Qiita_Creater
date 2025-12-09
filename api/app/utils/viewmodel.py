"""ViewModel変換ユーティリティ - ArticleState → ArticleViewModel"""
from typing import Dict, Any, Optional
from datetime import datetime


def phase_to_status_text(phase: str) -> str:
    """Phase名を自然文に変換"""
    phase_map = {
        "plan": "調査・構成設計中",
        "do": "執筆・検証・レビュー中",
        "check": "パフォーマンス分析中",
        "act": "改善指針の抽象化中",
        "completed": "完了"
    }
    return phase_map.get(phase, phase)


def get_next_action_hint(phase: str, article_state: Dict[str, Any]) -> str:
    """次の予定アクションを生成"""
    if phase == "plan":
        return "次：記事本文の執筆"
    elif phase == "do":
        review_result = article_state.get("review_result", {})
        if review_result.get("approval"):
            return "次：Qiitaへの投稿"
        else:
            return "次：レビュー結果に基づく修正"
    elif phase == "check":
        return "次：改善指針の抽象化"
    elif phase == "act":
        return "次：次の記事の計画"
    elif phase == "completed":
        return "完了済み"
    else:
        return "待機中"


def article_state_to_viewmodel(article_state: Dict[str, Any]) -> Dict[str, Any]:
    """ArticleStateをArticleViewModelに変換"""
    phase = article_state.get("current_phase", "plan")
    plan = article_state.get("plan", {})
    
    # titleの取得（planから、またはtopicから、またはデフォルト値）
    title = None
    if isinstance(plan, dict):
        title = plan.get("title")
    if not title:
        title = article_state.get("topic")
    if not title:
        title = "タイトル未設定"
    
    # markdownの取得（contentから、またはデフォルト値）
    markdown = article_state.get("content") or ""
    
    return {
        "id": article_state.get("article_id"),
        "title": title,
        "uiStatusText": phase_to_status_text(phase),
        "phase": phase,
        "markdown": markdown,
        "plan": plan,
        "reviewResult": article_state.get("review_result"),
        "qiitaUrl": article_state.get("qiita_url"),
        "qiitaItemId": article_state.get("qiita_item_id"),
        "kpiSummary": article_state.get("kpi"),
        "analysisResults": article_state.get("analysis_results"),
        "nextActionHint": get_next_action_hint(phase, article_state),
        "createdAt": article_state.get("created_at"),
        "updatedAt": article_state.get("updated_at"),
        "feedbackHistory": article_state.get("human_feedback", [])
    }

