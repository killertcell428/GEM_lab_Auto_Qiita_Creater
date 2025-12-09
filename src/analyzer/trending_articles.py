"""トレンド記事取得モジュール - Qiita APIで伸びている記事を取得"""
import requests
from typing import Dict, Any, List, Optional
from datetime import datetime, timedelta
from pathlib import Path
import json
from src.config_loader import get_config
from src.publisher.qiita_fetcher import QiitaFetcher


class TrendingArticlesFetcher:
    """Qiitaで伸びている記事を取得するクラス"""
    
    def __init__(self):
        """初期化"""
        config = get_config()
        qiita_config = config.get("api", {}).get("qiita", {})
        self.base_url = qiita_config.get("base_url", "https://qiita.com/api/v2")
        self.access_token = qiita_config.get("access_token", "")
        self.timeout = qiita_config.get("timeout", 30)
        self.qiita_fetcher = QiitaFetcher()
        
        # データ保存ディレクトリ
        self.data_dir = Path("data/analysis")
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.trending_file = self.data_dir / "trending_articles.json"
    
    def fetch_trending_by_likes(self, limit: int = 20, days: int = 7) -> List[Dict[str, Any]]:
        """
        いいね数が多い記事を取得
        
        Args:
            limit: 取得する記事数
            days: 過去何日分の記事を対象にするか
            
        Returns:
            List[Dict[str, Any]]: トレンド記事のリスト
        """
        headers = {
            "Authorization": f"Bearer {self.access_token}",
            "Content-Type": "application/json"
        } if self.access_token else {}
        
        cutoff_date = datetime.now() - timedelta(days=days)
        
        try:
            # Qiita APIで記事を取得（ページネーション対応）
            all_items = []
            page = 1
            per_page = 100
            
            while len(all_items) < limit:
                params = {
                    "page": page,
                    "per_page": min(per_page, limit - len(all_items))
                }
                
                response = requests.get(
                    f"{self.base_url}/items",
                    headers=headers,
                    params=params,
                    timeout=self.timeout
                )
                response.raise_for_status()
                items = response.json()
                
                if not items:
                    break
                
                # 日付でフィルタリング
                for item in items:
                    created_at = datetime.fromisoformat(
                        item.get("created_at", "").replace("Z", "+00:00")
                    )
                    if created_at.replace(tzinfo=None) >= cutoff_date:
                        all_items.append(item)
                
                if len(items) < per_page:
                    break
                
                page += 1
                
                # レート制限を考慮して少し待機
                import time
                time.sleep(0.5)
            
            # いいね数でソート
            all_items.sort(key=lambda x: x.get("likes_count", 0), reverse=True)
            
            return all_items[:limit]
            
        except requests.exceptions.RequestException as e:
            print(f"[WARN] トレンド記事の取得に失敗: {str(e)}")
            return []
    
    def fetch_trending_by_tag(self, tag: str, limit: int = 20, days: int = 7) -> List[Dict[str, Any]]:
        """
        特定タグの人気記事を取得
        
        Args:
            tag: タグ名
            limit: 取得する記事数
            days: 過去何日分の記事を対象にするか
            
        Returns:
            List[Dict[str, Any]]: トレンド記事のリスト
        """
        headers = {
            "Authorization": f"Bearer {self.access_token}",
            "Content-Type": "application/json"
        } if self.access_token else {}
        
        cutoff_date = datetime.now() - timedelta(days=days)
        
        try:
            all_items = []
            page = 1
            per_page = 100
            
            while len(all_items) < limit:
                params = {
                    "page": page,
                    "per_page": min(per_page, limit - len(all_items)),
                    "query": f"tag:{tag}"
                }
                
                response = requests.get(
                    f"{self.base_url}/items",
                    headers=headers,
                    params=params,
                    timeout=self.timeout
                )
                response.raise_for_status()
                items = response.json()
                
                if not items:
                    break
                
                # 日付でフィルタリング
                for item in items:
                    created_at = datetime.fromisoformat(
                        item.get("created_at", "").replace("Z", "+00:00")
                    )
                    if created_at.replace(tzinfo=None) >= cutoff_date:
                        all_items.append(item)
                
                if len(items) < per_page:
                    break
                
                page += 1
                
                import time
                time.sleep(0.5)
            
            # いいね数でソート
            all_items.sort(key=lambda x: x.get("likes_count", 0), reverse=True)
            
            return all_items[:limit]
            
        except requests.exceptions.RequestException as e:
            print(f"[WARN] タグ '{tag}' のトレンド記事取得に失敗: {str(e)}")
            return []
    
    def fetch_my_popular_articles(self, limit: int = 10, days: int = 90) -> List[Dict[str, Any]]:
        """
        自分の投稿から閲覧数・いいね数が高い記事を取得
        
        Args:
            limit: 取得する記事数
            days: 過去何日分の記事を対象にするか
            
        Returns:
            List[Dict[str, Any]]: 人気記事のリスト
        """
        try:
            items = self.qiita_fetcher.fetch_all_my_items(days_back=days)
            
            # いいね数と閲覧数の合計でソート
            items.sort(
                key=lambda x: x.get("likes_count", 0) + x.get("page_views_count", 0) * 0.1,
                reverse=True
            )
            
            return items[:limit]
            
        except Exception as e:
            print(f"[WARN] 自分の人気記事の取得に失敗: {str(e)}")
            return []
    
    def fetch_all_trending(self, tags: Optional[List[str]] = None, limit_per_source: int = 20) -> List[Dict[str, Any]]:
        """
        全てのソースからトレンド記事を取得
        
        Args:
            tags: 取得するタグのリスト（Noneの場合はデフォルトタグを使用）
            limit_per_source: 各ソースから取得する記事数
            
        Returns:
            List[Dict[str, Any]]: トレンド記事のリスト（重複除去済み）
        """
        all_articles = []
        article_ids = set()
        
        # 1. いいね数が多い記事
        print("[INFO] いいね数が多い記事を取得中...")
        trending_by_likes = self.fetch_trending_by_likes(limit=limit_per_source, days=7)
        for article in trending_by_likes:
            article_id = article.get("id")
            if article_id and article_id not in article_ids:
                article_ids.add(article_id)
                article["_source"] = "trending_by_likes"
                all_articles.append(article)
        
        # 2. 特定タグの人気記事
        if tags:
            for tag in tags:
                print(f"[INFO] タグ '{tag}' の人気記事を取得中...")
                trending_by_tag = self.fetch_trending_by_tag(tag, limit=limit_per_source, days=7)
                for article in trending_by_tag:
                    article_id = article.get("id")
                    if article_id and article_id not in article_ids:
                        article_ids.add(article_id)
                        article["_source"] = f"trending_by_tag_{tag}"
                        all_articles.append(article)
        
        # 3. 自分の人気記事
        print("[INFO] 自分の人気記事を取得中...")
        my_popular = self.fetch_my_popular_articles(limit=limit_per_source, days=90)
        for article in my_popular:
            article_id = article.get("id")
            if article_id and article_id not in article_ids:
                article_ids.add(article_id)
                article["_source"] = "my_popular"
                all_articles.append(article)
        
        # データを保存
        self.save_trending_articles(all_articles)
        
        return all_articles
    
    def save_trending_articles(self, articles: List[Dict[str, Any]]) -> None:
        """
        トレンド記事をJSONファイルに保存
        
        Args:
            articles: 保存する記事のリスト
        """
        data = {
            "fetched_at": datetime.now().isoformat(),
            "articles": articles
        }
        
        try:
            with open(self.trending_file, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, indent=2)
            print(f"[OK] {len(articles)}件のトレンド記事を保存しました: {self.trending_file}")
        except Exception as e:
            print(f"[WARN] トレンド記事の保存に失敗: {str(e)}")
    
    def load_trending_articles(self) -> List[Dict[str, Any]]:
        """
        保存されたトレンド記事を読み込む
        
        Returns:
            List[Dict[str, Any]]: トレンド記事のリスト
        """
        if not self.trending_file.exists():
            return []
        
        try:
            with open(self.trending_file, "r", encoding="utf-8") as f:
                data = json.load(f)
            return data.get("articles", [])
        except Exception as e:
            print(f"[WARN] トレンド記事の読み込みに失敗: {str(e)}")
            return []

