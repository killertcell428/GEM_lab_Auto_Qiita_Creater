"""いいね数の多い記事を分析して特徴を抽出するスクリプト"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.analyzer.trending_articles import TrendingArticlesFetcher
from src.analyzer.article_analyzer import ArticleAnalyzer
import json

def main():
    print("いいね数の多い記事を分析中...")
    
    # トレンド記事を取得
    fetcher = TrendingArticlesFetcher()
    articles = fetcher.fetch_trending_by_likes(limit=20, days=30)
    
    print(f"取得した記事数: {len(articles)}")
    
    if not articles:
        print("記事が取得できませんでした")
        return
    
    # 記事を分析
    analyzer = ArticleAnalyzer()
    analyzed = analyzer.analyze_multiple_articles(articles)
    
    # 共通特徴を抽出
    features = analyzer.extract_common_features(analyzed, top_n=10)
    
    # 保存
    analyzer.save_features(features)
    
    print("\n=== 分析結果 ===")
    print(json.dumps(features, ensure_ascii=False, indent=2))
    
    # 上位記事のタイトルと特徴を表示
    print("\n=== 上位記事のタイトル ===")
    sorted_articles = sorted(
        analyzed,
        key=lambda x: x.get("likes_count", 0) + x.get("page_views_count", 0) * 0.1,
        reverse=True
    )
    for i, article in enumerate(sorted_articles[:5], 1):
        print(f"\n{i}. {article.get('article_id', 'N/A')}")
        print(f"   いいね数: {article.get('likes_count', 0)}")
        print(f"   閲覧数: {article.get('page_views_count', 0)}")
        print(f"   構造: {article.get('structure', 'N/A')[:100]}...")
        print(f"   重要なポイント: {article.get('key_points', [])[:3]}")

if __name__ == "__main__":
    main()

