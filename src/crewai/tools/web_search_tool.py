"""Web検索ツール - CrewAI Toolとして実装"""
try:
    from crewai_tools import tool
except ImportError:
    from crewai.tools import tool
from typing import Dict, Any, List, Optional
import requests
from src.config_loader import get_config


@tool("Web検索ツール")
def search_web(
    query: str,
    num_results: int = 5,
    search_type: str = "general"
) -> str:
    """
    Webで検索を行い、公式ドキュメントやブログ記事を取得する
    
    Args:
        query: 検索クエリ（例: "Python pandas 公式ドキュメント"）
        num_results: 取得する結果数（デフォルト: 5）
        search_type: 検索タイプ（"general", "documentation", "blog", "github"）
    
    Returns:
        str: 検索結果の要約（JSON形式の文字列）
    """
    config = get_config()
    
    # SerperAPIまたはTavily APIを使用（設定ファイルから取得）
    # デフォルトではGemini APIの検索機能を使用
    try:
        # まず、Gemini APIの検索機能を試す
        from src.crewai.tools.gemini_tool import call_gemini_api
        
        search_prompt = f"""
以下のクエリでWeb検索を行い、公式ドキュメントやブログ記事を探してください：
検索クエリ: {query}
検索タイプ: {search_type}
取得数: {num_results}

検索結果を以下のJSON形式で返してください：
{{
    "query": "{query}",
    "results": [
        {{
            "title": "記事タイトル",
            "url": "URL",
            "snippet": "記事の抜粋",
            "type": "documentation|blog|github|other"
        }}
    ],
    "summary": "検索結果の要約"
}}
"""
        
        result = call_gemini_api(search_prompt)
        return result
        
    except Exception as e:
        # フォールバック: 検索結果の代わりにエラーメッセージを返す
        return f'{{"error": "Web検索に失敗しました: {str(e)}", "query": "{query}"}}'


@tool("公式ドキュメント検索ツール")
def search_official_documentation(
    library_name: str,
    topic: str = ""
) -> str:
    """
    特定のライブラリの公式ドキュメントを検索する
    
    Args:
        library_name: ライブラリ名（例: "pandas", "numpy", "biopython"）
        topic: 検索トピック（例: "CNV解析", "データ読み込み"）
    
    Returns:
        str: 公式ドキュメントの情報（JSON形式の文字列）
    """
    query = f"{library_name} 公式ドキュメント {topic}".strip()
    return search_web(query, num_results=3, search_type="documentation")


@tool("実装例検索ツール")
def search_implementation_examples(
    technology: str,
    use_case: str = ""
) -> str:
    """
    実装例やブログ記事を検索する
    
    Args:
        technology: 技術名（例: "CNV解析", "トランスクリプトーム解析"）
        use_case: ユースケース（例: "Python実装", "バイオインフォマティクス"）
    
    Returns:
        str: 実装例の情報（JSON形式の文字列）
    """
    query = f"{technology} {use_case} 実装例 ブログ".strip()
    return search_web(query, num_results=5, search_type="blog")

