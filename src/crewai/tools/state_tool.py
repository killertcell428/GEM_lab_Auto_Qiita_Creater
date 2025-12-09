"""State管理ツール - CrewAI Toolとして実装"""
from crewai.tools import tool
from typing import Dict, Any, Optional
from src.crewai.state.article_state import ArticleState
from src.crewai.state.pdca_state import PDCAState


@tool
def load_article_state(article_id: str) -> str:
    """
    記事のStateを読み込む
    
    Args:
        article_id: 記事ID
    
    Returns:
        Article State（JSON形式の文字列）
    """
    try:
        state = ArticleState.load(article_id)
        import json
        return json.dumps(state.to_dict(), ensure_ascii=False)
    except Exception as e:
        raise Exception(f"Article Stateの読み込みに失敗しました: {str(e)}")


@tool
def save_article_state(article_id: str, state_data: Dict[str, Any]) -> str:
    """
    記事のStateを保存する
    
    Args:
        article_id: 記事ID
        state_data: Stateデータ（辞書形式）
    
    Returns:
        保存結果（JSON形式の文字列）
    """
    try:
        state = ArticleState(**{**state_data, "article_id": article_id})
        state.save()
        return f"Article Stateを保存しました: {article_id}"
    except Exception as e:
        raise Exception(f"Article Stateの保存に失敗しました: {str(e)}")


@tool
def load_pdca_state(cycle_id: str) -> Optional[str]:
    """
    PDCA Stateを読み込む
    
    Args:
        cycle_id: サイクルID
    
    Returns:
        PDCA State（JSON形式の文字列）、存在しない場合はNone
    """
    try:
        state = PDCAState.load(cycle_id)
        if state is None:
            return None
        import json
        return json.dumps(state.to_dict(), ensure_ascii=False)
    except Exception as e:
        raise Exception(f"PDCA Stateの読み込みに失敗しました: {str(e)}")


@tool
def save_pdca_state(cycle_id: str, state_data: Dict[str, Any]) -> str:
    """
    PDCA Stateを保存する
    
    Args:
        cycle_id: サイクルID
        state_data: Stateデータ（辞書形式）
    
    Returns:
        保存結果（JSON形式の文字列）
    """
    try:
        state = PDCAState(**{**state_data, "cycle_id": cycle_id})
        state.save()
        return f"PDCA Stateを保存しました: {cycle_id}"
    except Exception as e:
        raise Exception(f"PDCA Stateの保存に失敗しました: {str(e)}")

