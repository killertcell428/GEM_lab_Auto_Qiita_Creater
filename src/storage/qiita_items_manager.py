"""Qiita投稿管理モジュール - 過去の投稿を保存・管理"""
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional


class QiitaItemsManager:
    """Qiitaから取得した過去の投稿を管理するクラス"""
    
    def __init__(self, data_dir: str = "data"):
        """
        Args:
            data_dir: データを保存するディレクトリ
        """
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(exist_ok=True)
        self.items_file = self.data_dir / "qiita_items.json"
        self._items_cache: Optional[List[Dict[str, Any]]] = None
    
    def save_items(self, items: List[Dict[str, Any]]) -> None:
        """
        Qiitaから取得した投稿を保存
        
        Args:
            items: 投稿のリスト
        """
        data = {
            "last_updated": datetime.now().isoformat(),
            "items": items,
            "total_count": len(items)
        }
        
        with open(self.items_file, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        
        # キャッシュを更新
        self._items_cache = items
    
    def load_items(self) -> List[Dict[str, Any]]:
        """
        保存された投稿を読み込む
        
        Returns:
            List[Dict[str, Any]]: 投稿のリスト
        """
        if self._items_cache is not None:
            return self._items_cache
        
        if not self.items_file.exists():
            return []
        
        with open(self.items_file, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        items = data.get("items", [])
        self._items_cache = items
        return items
    
    def get_items_summary(self, limit: int = 10) -> str:
        """
        保存された投稿の要約を取得（記事生成プロンプト用）
        
        Args:
            limit: 取得する投稿数
            
        Returns:
            str: 投稿の要約テキスト
        """
        items = self.load_items()
        
        if not items:
            return "過去の投稿はありません。"
        
        # 最新の投稿を制限数まで取得
        recent_items = items[:limit]
        
        summary_lines = ["## 過去の投稿（参考）\n"]
        
        for i, item in enumerate(recent_items, 1):
            title = item.get("title", "タイトルなし")
            url = item.get("url", "")
            created_at = item.get("created_at", "")
            tags = [tag.get("name", "") if isinstance(tag, dict) else tag for tag in item.get("tags", [])]
            
            # 日付をフォーマット
            try:
                date_obj = datetime.fromisoformat(created_at.replace("Z", "+00:00"))
                date_str = date_obj.strftime("%Y年%m月%d日")
            except:
                date_str = created_at[:10] if created_at else "日付不明"
            
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
        
        all_items = self.load_items()
        
        # タグの一致度でソート
        scored_items = []
        for item in all_items:
            item_tags = []
            for tag in item.get("tags", []):
                if isinstance(tag, dict):
                    item_tags.append(tag.get("name", "").lower())
                else:
                    item_tags.append(str(tag).lower())
            
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
    
    def clear_cache(self) -> None:
        """キャッシュをクリア"""
        self._items_cache = None

