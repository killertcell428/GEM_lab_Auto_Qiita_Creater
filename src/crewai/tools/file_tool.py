"""ファイルI/Oツール - CrewAI Toolとして実装"""
from crewai.tools import tool
from typing import Dict, Any, Optional
from pathlib import Path
from src.generator.draft_loader import load_draft_from_file, parse_draft


@tool
def load_draft_file(file_path: str) -> str:
    """
    ドラフトファイルを読み込む
    
    Args:
        file_path: ドラフトファイルのパス
    
    Returns:
        ドラフト内容（Markdown形式）
    """
    try:
        draft_data = load_draft_from_file(file_path)
        return draft_data.get("raw_content", draft_data.get("content", ""))
    except Exception as e:
        raise Exception(f"ドラフトファイルの読み込みに失敗しました: {str(e)}")


@tool
def save_article_file(
    file_path: str,
    content: str,
    title: Optional[str] = None
) -> str:
    """
    記事をMarkdownファイルとして保存する
    
    Args:
        file_path: 保存先のファイルパス
        content: 記事内容（Markdown形式）
        title: 記事タイトル（オプション、ファイル名に使用）
    
    Returns:
        保存結果メッセージ
    """
    try:
        file_path_obj = Path(file_path)
        file_path_obj.parent.mkdir(parents=True, exist_ok=True)
        
        # タイトルがある場合は先頭に追加
        if title:
            content = f"# {title}\n\n{content}"
        
        with open(file_path_obj, "w", encoding="utf-8") as f:
            f.write(content)
        
        return f"記事を保存しました: {file_path}"
    except Exception as e:
        raise Exception(f"記事ファイルの保存に失敗しました: {str(e)}")


@tool
def load_article_file(file_path: str) -> str:
    """
    記事ファイルを読み込む
    
    Args:
        file_path: 記事ファイルのパス
    
    Returns:
        記事内容（Markdown形式）
    """
    try:
        file_path_obj = Path(file_path)
        if not file_path_obj.exists():
            raise FileNotFoundError(f"ファイルが見つかりません: {file_path}")
        
        with open(file_path_obj, "r", encoding="utf-8") as f:
            content = f.read()
        
        return content
    except Exception as e:
        raise Exception(f"記事ファイルの読み込みに失敗しました: {str(e)}")

