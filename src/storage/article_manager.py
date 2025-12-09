"""記事管理モジュール - 記事の保存・管理"""
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List, Optional
from enum import Enum


class ArticleStatus(Enum):
    """記事のステータス"""
    DRAFT = "draft"
    APPROVED = "approved"
    PUBLISHED = "published"


class ArticleManager:
    """記事の保存・管理を行うクラス"""
    
    def __init__(self, articles_dir: str = "articles", data_dir: str = "data"):
        """
        Args:
            articles_dir: 記事ファイル（Markdown）を保存するディレクトリ
            data_dir: メタデータ（JSON）を保存するディレクトリ
        """
        self.articles_dir = Path(articles_dir)
        self.data_dir = Path(data_dir)
        
        # ディレクトリが存在しない場合は作成
        self.articles_dir.mkdir(exist_ok=True)
        self.data_dir.mkdir(exist_ok=True)
    
    def save_article(self, article_data: Dict[str, Any], status: ArticleStatus = ArticleStatus.DRAFT) -> str:
        """
        記事を保存
        
        Args:
            article_data: 記事データ（title, content, tagsを含む）
            status: 記事のステータス
            
        Returns:
            str: 記事ID（ファイル名のベース）
        """
        # 記事IDを生成（タイムスタンプベース）
        article_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Markdownファイルとして保存
        title = article_data.get("title", "untitled")
        # ファイル名に使用できない文字を置換
        safe_title = "".join(c if c.isalnum() or c in (" ", "-", "_") else "_" for c in title)
        safe_title = safe_title[:50]  # 長さ制限
        
        md_filename = f"{article_id}_{safe_title}.md"
        md_path = self.articles_dir / md_filename
        
        # Markdownコンテンツを構築
        md_content = f"""# {title}

{article_data.get('content', '')}

---
生成日時: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
ステータス: {status.value}
タグ: {', '.join(article_data.get('tags', []))}
"""
        
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(md_content)
        
        # メタデータをJSONで保存
        metadata = {
            "article_id": article_id,
            "title": title,
            "status": status.value,
            "tags": article_data.get("tags", []),
            "created_at": datetime.now().isoformat(),
            "updated_at": datetime.now().isoformat(),
            "md_filename": md_filename,
            "scheduled_publish_at": None,
            "published_at": None,
            "qiita_item_id": None
        }
        
        metadata_filename = f"{article_id}_metadata.json"
        metadata_path = self.data_dir / metadata_filename
        
        with open(metadata_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, ensure_ascii=False, indent=2)
        
        return article_id
    
    def load_article(self, article_id: str) -> Dict[str, Any]:
        """
        記事を読み込む
        
        Args:
            article_id: 記事ID
            
        Returns:
            Dict[str, Any]: 記事データとメタデータ
        """
        # メタデータを読み込む
        metadata_pattern = f"{article_id}_metadata.json"
        metadata_files = list(self.data_dir.glob(metadata_pattern))
        
        if not metadata_files:
            raise FileNotFoundError(f"記事ID {article_id} が見つかりません")
        
        with open(metadata_files[0], "r", encoding="utf-8") as f:
            metadata = json.load(f)
        
        # Markdownファイルを読み込む
        md_filename = metadata.get("md_filename")
        if not md_filename:
            raise ValueError(f"メタデータにmd_filenameがありません: {article_id}")
        
        md_path = self.articles_dir / md_filename
        if not md_path.exists():
            raise FileNotFoundError(f"Markdownファイルが見つかりません: {md_filename}")
        
        with open(md_path, "r", encoding="utf-8") as f:
            content = f.read()
        
        # タイトルとコンテンツを抽出（簡単なパース）
        lines = content.split("\n")
        title = ""
        body_lines = []
        found_title = False
        
        for i, line in enumerate(lines):
            if line.startswith("# ") and not found_title:
                title = line.replace("# ", "").strip()
                found_title = True
            elif line.startswith("---"):
                # メタデータセクションの開始（スキップ）
                continue
            elif found_title:
                # タイトルが見つかった後の行は本文
                body_lines.append(line)
        
        content_body = "\n".join(body_lines).strip()
        
        # タイトルが取得できなかった場合、メタデータから取得
        if not title:
            title = metadata.get("title", "")
        
        return {
            "article_id": article_id,
            "title": title or metadata.get("title", ""),
            "content": content_body,
            "tags": metadata.get("tags", []),
            "status": metadata.get("status", "draft"),
            "metadata": metadata
        }
    
    def update_status(self, article_id: str, status: ArticleStatus) -> None:
        """
        記事のステータスを更新
        
        Args:
            article_id: 記事ID
            status: 新しいステータス
        """
        metadata_pattern = f"{article_id}_metadata.json"
        metadata_files = list(self.data_dir.glob(metadata_pattern))
        
        if not metadata_files:
            raise FileNotFoundError(f"記事ID {article_id} が見つかりません")
        
        with open(metadata_files[0], "r", encoding="utf-8") as f:
            metadata = json.load(f)
        
        metadata["status"] = status.value
        metadata["updated_at"] = datetime.now().isoformat()
        
        if status == ArticleStatus.PUBLISHED:
            metadata["published_at"] = datetime.now().isoformat()
        
        with open(metadata_files[0], "w", encoding="utf-8") as f:
            json.dump(metadata, f, ensure_ascii=False, indent=2)
    
    def list_articles(self, status: Optional[ArticleStatus] = None) -> List[Dict[str, Any]]:
        """
        記事のリストを取得
        
        Args:
            status: フィルタするステータス（Noneの場合は全て）
            
        Returns:
            List[Dict[str, Any]]: 記事のメタデータのリスト
        """
        metadata_files = list(self.data_dir.glob("*_metadata.json"))
        articles = []
        
        for metadata_file in metadata_files:
            with open(metadata_file, "r", encoding="utf-8") as f:
                metadata = json.load(f)
            
            if status is None or metadata.get("status") == status.value:
                articles.append(metadata)
        
        # 作成日時の降順でソート
        articles.sort(key=lambda x: x.get("created_at", ""), reverse=True)
        
        return articles
    
    def get_approved_articles(self) -> List[Dict[str, Any]]:
        """承認済み記事のリストを取得"""
        return self.list_articles(ArticleStatus.APPROVED)

