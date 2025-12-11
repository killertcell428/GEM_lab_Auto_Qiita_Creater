"""ViewModel変換ユーティリティ - ArticleState → ArticleViewModel"""
from typing import Dict, Any, Optional
from datetime import datetime
import json


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


def is_json_string(value: str) -> bool:
    """文字列がJSON形式かどうかを判定"""
    if not isinstance(value, str) or not value.strip():
        return False
    try:
        json.loads(value)
        return True
    except (json.JSONDecodeError, ValueError):
        return False


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
    content = article_state.get("content") or ""
    
    # JSON文字列が入っている場合は警告を出す
    if is_json_string(content):
        article_id = article_state.get("article_id", "unknown")
        print(f"[WARN] article_id={article_id}: contentフィールドにJSONが入っています。Markdownコンテンツが表示されません。")
        # 空文字列にフォールバック（ユーザーにエラーメッセージを表示するため）
        markdown = ""
    else:
        markdown = content
    
    # 承認待ちの場合のヒントを更新
    next_hint = get_next_action_hint(phase, article_state)
    if article_state.get("pending_approval") and article_state.get("approval_status") == "pending":
        deadline = article_state.get("approval_deadline")
        if deadline:
            try:
                deadline_dt = datetime.fromisoformat(deadline)
                now = datetime.now()
                if now < deadline_dt:
                    remaining = deadline_dt - now
                    hours = int(remaining.total_seconds() / 3600)
                    minutes = int((remaining.total_seconds() % 3600) / 60)
                    next_hint = f"承認待ち（残り時間: {hours}時間{minutes}分）"
                else:
                    next_hint = "承認期限超過（自動投稿予定）"
            except:
                next_hint = "承認待ち"
    
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
        "nextActionHint": next_hint,
        "createdAt": article_state.get("created_at"),
        "updatedAt": article_state.get("updated_at"),
        "feedbackHistory": article_state.get("human_feedback", []),
        # 承認関連フィールド
        "pendingApproval": article_state.get("pending_approval", False),
        "approvalDeadline": article_state.get("approval_deadline"),
        "approvalStatus": article_state.get("approval_status"),
        "scheduledPublishDate": article_state.get("scheduled_publish_date"),
        # リサーチ結果
        "researchReport": article_state.get("research_report")
    }

