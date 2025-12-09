"""Qiita API呼び出しツール - CrewAI Toolとして実装"""
from crewai.tools import tool
from typing import Dict, Any, List, Optional
from src.publisher.qiita_fetcher import QiitaFetcher
from src.publisher.qiita_publisher import QiitaPublisher


@tool
def fetch_qiita_articles(
    days_back: int = 180,
    limit: Optional[int] = None
) -> List[Dict[str, Any]]:
    """
    Qiitaから過去の記事を取得する
    
    Args:
        days_back: 過去何日分の記事を取得するか（デフォルト: 180）
        limit: 取得する記事数の上限（Noneの場合は全て）
    
    Returns:
        記事のリスト（JSON形式の文字列）
    """
    try:
        fetcher = QiitaFetcher()
        items = fetcher.fetch_all_my_items(days_back=days_back)
        
        if limit:
            items = items[:limit]
        
        # JSON形式の文字列として返す
        import json
        return json.dumps(items, ensure_ascii=False)
    except Exception as e:
        raise Exception(f"Qiita記事の取得に失敗しました: {str(e)}")


@tool
def publish_to_qiita(
    title: str,
    content: str,
    tags: List[str],
    private: bool = False,
    tweet: bool = False
) -> Dict[str, Any]:
    """
    Qiitaに記事を投稿する
    
    Args:
        title: 記事タイトル
        content: 記事本文（Markdown形式）
        tags: タグのリスト（最大5個）
        private: 限定公開にするかどうか
        tweet: ツイートするかどうか
    
    Returns:
        投稿結果（JSON形式の文字列）
    """
    try:
        publisher = QiitaPublisher()
        
        article_data = {
            "title": title,
            "content": content,
            "tags": tags[:5]  # Qiita APIの制限
        }
        
        result = publisher.publish_article(
            article_data,
            private=private,
            tweet=tweet
        )
        
        # JSON形式の文字列として返す
        import json
        return json.dumps(result, ensure_ascii=False)
    except Exception as e:
        raise Exception(f"Qiitaへの投稿に失敗しました: {str(e)}")


@tool
def fetch_trending_articles(
    limit: int = 20,
    days: int = 7
) -> str:
    """
    Qiitaでトレンドの記事を取得する
    
    Args:
        limit: 取得する記事数（デフォルト: 20）
        days: 過去何日分の記事を対象にするか（デフォルト: 7）
    
    Returns:
        トレンド記事のリスト（JSON形式の文字列）
    """
    try:
        from src.analyzer.trending_articles import TrendingArticlesFetcher
        
        fetcher = TrendingArticlesFetcher()
        articles = fetcher.fetch_all_trending(
            tags=None,
            limit_per_source=limit
        )
        
        # 制限を適用
        if limit:
            articles = articles[:limit]
        
        # JSON形式の文字列として返す
        import json
        return json.dumps(articles, ensure_ascii=False)
    except Exception as e:
        raise Exception(f"トレンド記事の取得に失敗しました: {str(e)}")

