"""フィードバック収集モジュール - 投稿後の記事の閲覧数・いいね数を取得"""
from typing import Dict, Any, List, Optional
from pathlib import Path
import json
from datetime import datetime
from src.publisher.qiita_fetcher import QiitaFetcher


class FeedbackCollector:
    """投稿後の記事のフィードバックを収集するクラス"""
    
    def __init__(self):
        """初期化"""
        self.qiita_fetcher = QiitaFetcher()
        
        # データ保存ディレクトリ
        self.data_dir = Path("data/analysis")
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.feedback_file = self.data_dir / "feedback_history.json"
    
    def collect_feedback_for_article(self, article_id: str) -> Optional[Dict[str, Any]]:
        """
        特定の記事のフィードバックを収集
        
        Args:
            article_id: Qiita記事ID
            
        Returns:
            Optional[Dict[str, Any]]: フィードバックデータ（取得できない場合はNone）
        """
        try:
            # 自分の投稿から該当記事を検索
            items = self.qiita_fetcher.fetch_all_my_items(days_back=365)
            
            for item in items:
                if item.get("id") == article_id:
                    return {
                        "article_id": article_id,
                        "title": item.get("title", ""),
                        "url": item.get("url", ""),
                        "likes_count": item.get("likes_count", 0),
                        "page_views_count": item.get("page_views_count", 0),
                        "comments_count": item.get("comments_count", 0),
                        "stocks_count": item.get("stocks_count", 0),
                        "created_at": item.get("created_at", ""),
                        "updated_at": item.get("updated_at", ""),
                        "collected_at": datetime.now().isoformat()
                    }
            
            return None
            
        except Exception as e:
            print(f"[WARN] フィードバック収集に失敗: {str(e)}")
            return None
    
    def collect_feedback_for_recent_articles(self, days: int = 30) -> List[Dict[str, Any]]:
        """
        最近投稿した記事のフィードバックを一括収集
        
        Args:
            days: 過去何日分の記事を対象にするか
            
        Returns:
            List[Dict[str, Any]]: フィードバックデータのリスト
        """
        try:
            items = self.qiita_fetcher.fetch_all_my_items(days_back=days)
            
            feedback_list = []
            for item in items:
                feedback = {
                    "article_id": item.get("id", ""),
                    "title": item.get("title", ""),
                    "url": item.get("url", ""),
                    "likes_count": item.get("likes_count", 0),
                    "page_views_count": item.get("page_views_count", 0),
                    "comments_count": item.get("comments_count", 0),
                    "stocks_count": item.get("stocks_count", 0),
                    "created_at": item.get("created_at", ""),
                    "updated_at": item.get("updated_at", ""),
                    "collected_at": datetime.now().isoformat()
                }
                feedback_list.append(feedback)
            
            return feedback_list
            
        except Exception as e:
            print(f"[WARN] フィードバック一括収集に失敗: {str(e)}")
            return []
    
    def save_feedback(self, feedback: Dict[str, Any]) -> None:
        """
        フィードバックを履歴に追加
        
        Args:
            feedback: 保存するフィードバックデータ
        """
        # 既存の履歴を読み込む
        history = self.load_feedback_history()
        
        # 同じ記事IDの既存データを更新または追加
        article_id = feedback.get("article_id")
        if article_id:
            # 既存のエントリを探す
            existing_index = None
            for i, entry in enumerate(history):
                if entry.get("article_id") == article_id:
                    existing_index = i
                    break
            
            if existing_index is not None:
                # 既存のエントリを更新（より新しいデータで上書き）
                history[existing_index] = feedback
            else:
                # 新しいエントリを追加
                history.append(feedback)
        
        # 履歴を保存
        data = {
            "updated_at": datetime.now().isoformat(),
            "history": history
        }
        
        try:
            with open(self.feedback_file, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, indent=2)
            print(f"[OK] フィードバックを保存しました: {self.feedback_file}")
        except Exception as e:
            print(f"[WARN] フィードバックの保存に失敗: {str(e)}")
    
    def load_feedback_history(self) -> List[Dict[str, Any]]:
        """
        フィードバック履歴を読み込む
        
        Returns:
            List[Dict[str, Any]]: フィードバック履歴
        """
        if not self.feedback_file.exists():
            return []
        
        try:
            with open(self.feedback_file, "r", encoding="utf-8") as f:
                data = json.load(f)
            return data.get("history", [])
        except Exception as e:
            print(f"[WARN] フィードバック履歴の読み込みに失敗: {str(e)}")
            return []
    
    def get_high_performance_articles(self, min_likes: int = 5, min_views: int = 100) -> List[Dict[str, Any]]:
        """
        パフォーマンスが高い記事を取得
        
        Args:
            min_likes: 最小いいね数
            min_views: 最小閲覧数
            
        Returns:
            List[Dict[str, Any]]: 高パフォーマンス記事のリスト
        """
        history = self.load_feedback_history()
        
        high_performance = []
        for entry in history:
            likes = entry.get("likes_count", 0)
            views = entry.get("page_views_count", 0)
            
            if likes >= min_likes and views >= min_views:
                high_performance.append(entry)
        
        # いいね数と閲覧数の合計でソート
        high_performance.sort(
            key=lambda x: x.get("likes_count", 0) + x.get("page_views_count", 0) * 0.1,
            reverse=True
        )
        
        return high_performance

