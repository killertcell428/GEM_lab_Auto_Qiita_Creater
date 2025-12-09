"""Qiita投稿取得モジュール - Qiita API v2を使用して過去の投稿を取得"""
import requests
from typing import Dict, Any, List, Optional
from datetime import datetime, timedelta
from src.config_loader import get_config


class QiitaFetcher:
    """Qiitaから過去の投稿を取得するクラス"""
    
    def __init__(self):
        """初期化"""
        config = get_config()
        qiita_config = config.get("api", {}).get("qiita", {})
        self.base_url = qiita_config.get("base_url", "https://qiita.com/api/v2")
        self.access_token = qiita_config.get("access_token", "")
        self.timeout = qiita_config.get("timeout", 30)
        
        if not self.access_token:
            raise ValueError(
                "QIITA_ACCESS_TOKENが設定されていません。.envファイルを確認してください。\n"
                "Qiita API アクセストークンは https://qiita.com/settings/applications で取得できます。"
            )
    
    def get_authenticated_user(self) -> Dict[str, Any]:
        """
        認証済みユーザー情報を取得
        
        Returns:
            Dict[str, Any]: ユーザー情報
        """
        headers = {
            "Authorization": f"Bearer {self.access_token}",
            "Content-Type": "application/json"
        }
        
        try:
            response = requests.get(
                f"{self.base_url}/authenticated_user",
                headers=headers,
                timeout=self.timeout
            )
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            raise Exception(f"ユーザー情報の取得に失敗しました: {str(e)}")
    
    def fetch_my_items(
        self,
        page: int = 1,
        per_page: int = 100,
        days_back: Optional[int] = None
    ) -> List[Dict[str, Any]]:
        """
        自分の投稿を取得
        
        Args:
            page: ページ番号
            per_page: 1ページあたりの件数（最大100）
            days_back: 過去何日分の投稿を取得するか（Noneの場合は全て）
            
        Returns:
            List[Dict[str, Any]]: 投稿のリスト
        """
        headers = {
            "Authorization": f"Bearer {self.access_token}",
            "Content-Type": "application/json"
        }
        
        # クエリパラメータを構築
        params = {
            "page": page,
            "per_page": min(per_page, 100)  # 最大100件
        }
        
        try:
            response = requests.get(
                f"{self.base_url}/authenticated_user/items",
                headers=headers,
                params=params,
                timeout=self.timeout
            )
            response.raise_for_status()
            items = response.json()
            
            # 日付でフィルタリング
            if days_back:
                cutoff_date = datetime.now() - timedelta(days=days_back)
                filtered_items = []
                for item in items:
                    created_at = datetime.fromisoformat(item.get("created_at", "").replace("Z", "+00:00"))
                    if created_at.replace(tzinfo=None) >= cutoff_date:
                        filtered_items.append(item)
                return filtered_items
            
            return items
        except requests.exceptions.RequestException as e:
            raise Exception(f"投稿の取得に失敗しました: {str(e)}")
    
    def fetch_all_my_items(self, days_back: Optional[int] = None) -> List[Dict[str, Any]]:
        """
        自分の全ての投稿を取得（ページネーション対応）
        
        Args:
            days_back: 過去何日分の投稿を取得するか（Noneの場合は全て）
            
        Returns:
            List[Dict[str, Any]]: 投稿のリスト
        """
        all_items = []
        page = 1
        per_page = 100
        
        while True:
            items = self.fetch_my_items(page=page, per_page=per_page, days_back=days_back)
            if not items:
                break
            
            all_items.extend(items)
            
            # 100件未満の場合は最後のページ
            if len(items) < per_page:
                break
            
            page += 1
        
        # 作成日時の降順でソート
        all_items.sort(
            key=lambda x: x.get("created_at", ""),
            reverse=True
        )
        
        return all_items
    
    def get_recent_items_summary(self, limit: int = 10) -> str:
        """
        最近の投稿の要約を取得（記事生成プロンプト用）
        
        Args:
            limit: 取得する投稿数
            
        Returns:
            str: 投稿の要約テキスト
        """
        items = self.fetch_all_my_items(days_back=90)  # 過去90日分
        
        if not items:
            return "過去の投稿はありません。"
        
        # 最新の投稿を制限数まで取得
        recent_items = items[:limit]
        
        summary_lines = ["## 過去の投稿（参考）\n"]
        
        for i, item in enumerate(recent_items, 1):
            title = item.get("title", "タイトルなし")
            url = item.get("url", "")
            created_at = item.get("created_at", "")
            tags = [tag.get("name", "") for tag in item.get("tags", [])]
            
            # 日付をフォーマット
            try:
                date_obj = datetime.fromisoformat(created_at.replace("Z", "+00:00"))
                date_str = date_obj.strftime("%Y年%m月%d日")
            except:
                date_str = created_at[:10]
            
            summary_lines.append(
                f"{i}. **[{title}]({url})** ({date_str})\n"
                f"   - タグ: {', '.join(tags[:5])}\n"
            )
        
        return "\n".join(summary_lines)
    
    def get_related_items(self, current_tags: List[str], limit: int = 5) -> List[Dict[str, Any]]:
        """
        現在の記事のタグに関連する過去の投稿を取得
        
        Args:
            current_tags: 現在の記事のタグ
            limit: 取得する投稿数
            
        Returns:
            List[Dict[str, Any]]: 関連投稿のリスト
        """
        if not current_tags:
            return []
        
        all_items = self.fetch_all_my_items(days_back=180)  # 過去180日分
        
        # タグの一致度でソート
        scored_items = []
        for item in all_items:
            item_tags = [tag.get("name", "").lower() for tag in item.get("tags", [])]
            current_tags_lower = [tag.lower() for tag in current_tags]
            
            # 共通タグの数をスコアとする
            common_tags = set(item_tags) & set(current_tags_lower)
            score = len(common_tags)
            
            if score > 0:
                scored_items.append((score, item))
        
        # スコアの降順でソート
        scored_items.sort(key=lambda x: x[0], reverse=True)
        
        # 上位の投稿を返す
        return [item for _, item in scored_items[:limit]]

