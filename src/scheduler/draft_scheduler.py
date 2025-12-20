"""ドラフト記事のスケジュール投稿管理"""
from pathlib import Path
from typing import Dict, Any, Optional, List
from datetime import datetime
import json
import sys

# プロジェクトルートをパスに追加
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from src.generator.draft_loader import load_draft_from_file
from src.publisher.qiita_publisher import QiitaPublisher
from src.config_loader import get_config
from src.scheduler.draft_processor import DraftProcessor


class DraftScheduler:
    """ドラフト記事のスケジュール投稿を管理するクラス"""
    
    def __init__(self, drafts_dir: Optional[str] = None):
        """
        初期化
        
        Args:
            drafts_dir: ドラフトファイルのディレクトリパス（デフォルト: data/drafts/docs/blog）
        """
        if drafts_dir is None:
            self.drafts_dir = project_root / "data" / "drafts" / "docs" / "blog"
        else:
            self.drafts_dir = Path(drafts_dir)
        
        self.state_file = project_root / "data" / "state" / "scheduled_posts.json"
        self.publisher = QiitaPublisher()
        self.config = get_config()
        self.processor = DraftProcessor()
        
        # 状態ファイルのディレクトリを作成
        self.state_file.parent.mkdir(parents=True, exist_ok=True)
    
    def _load_state(self) -> Dict[str, Any]:
        """投稿状態を読み込む"""
        if not self.state_file.exists():
            return {
                "published_drafts": [],
                "last_published_date": None,
                "current_index": 0
            }
        
        try:
            with open(self.state_file, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            print(f"[WARN] 状態ファイルの読み込みエラー: {str(e)}")
            return {
                "published_drafts": [],
                "last_published_date": None,
                "current_index": 0
            }
    
    def _save_state(self, state: Dict[str, Any]):
        """投稿状態を保存"""
        try:
            with open(self.state_file, "w", encoding="utf-8") as f:
                json.dump(state, f, ensure_ascii=False, indent=2)
        except Exception as e:
            print(f"[ERROR] 状態ファイルの保存エラー: {str(e)}")
    
    def get_draft_files(self) -> List[str]:
        """ドラフトファイルのリストを取得（順番に）"""
        # 01_introduction.md から 08_summary.md まで順番に取得
        expected_files = [
            "01_introduction.md",
            "02_environment_setup.md",
            "03_simulation.md",
            "04_alignment_counting.md",
            "05_deseq2_analysis.md",
            "06_visualization.md",
            "07_cnv_analysis.md",
            "08_summary.md"
        ]
        
        draft_files = []
        for filename in expected_files:
            file_path = self.drafts_dir / filename
            if file_path.exists():
                draft_files.append(str(file_path))
            else:
                print(f"[WARN] ドラフトファイルが見つかりません: {filename}")
        
        return draft_files
    
    def get_next_draft(self) -> Optional[str]:
        """次に投稿すべきドラフトファイルを取得"""
        state = self._load_state()
        draft_files = self.get_draft_files()
        
        if not draft_files:
            print("[WARN] ドラフトファイルが見つかりません")
            return None
        
        current_index = state.get("current_index", 0)
        
        # すべて投稿済みの場合、最初から再開
        if current_index >= len(draft_files):
            print("[INFO] すべての記事を投稿済みです。最初から再開します。")
            current_index = 0
            state["current_index"] = 0
            state["published_drafts"] = []
            self._save_state(state)
        
        return draft_files[current_index]
    
    def publish_next_draft(self, auto_publish: bool = True) -> Dict[str, Any]:
        """
        次のドラフトを投稿
        
        Args:
            auto_publish: 自動的に投稿するかどうか（デフォルト: True）
            
        Returns:
            Dict[str, Any]: 投稿結果
        """
        try:
            # 次のドラフトファイルを取得
            draft_file = self.get_next_draft()
            if not draft_file:
                return {
                    "success": False,
                    "error": "投稿するドラフトファイルが見つかりません"
                }
            
            print(f"[SCHEDULER] ドラフトファイルを読み込み: {draft_file}")
            
            # ドラフトを読み込む
            draft_data = load_draft_from_file(draft_file)
            
            title = draft_data.get("title", "タイトル未設定")
            raw_content = draft_data.get("raw_content", draft_data.get("content", ""))
            
            # 記事内容を処理（パス修正、画像・テーブル処理）
            print(f"[SCHEDULER] 記事内容を処理中...")
            content = self.processor.process_content(raw_content, Path(draft_file))
            
            tags = draft_data.get("tags", [])
            
            # デフォルトタグを追加
            config_tags = self.config.get("article", {}).get("default_tags", [])
            if config_tags:
                # 既存のタグと重複しないように追加
                for tag in config_tags:
                    if tag not in tags:
                        tags.append(tag)
            
            # バイオインフォマティクス関連のタグを追加
            bio_tags = ["バイオインフォマティクス", "RNA-seq", "Python"]
            for tag in bio_tags:
                if tag not in tags:
                    tags.append(tag)
            
            # タグ数を制限（最大5個）
            tags = tags[:5]
            
            print(f"[SCHEDULER] タイトル: {title}")
            print(f"[SCHEDULER] タグ: {tags}")
            
            if not auto_publish:
                return {
                    "success": True,
                    "draft_file": draft_file,
                    "title": title,
                    "content": content[:200] + "..." if len(content) > 200 else content,
                    "tags": tags,
                    "published": False,
                    "message": "承認待ち（自動投稿は無効）"
                }
            
            # Qiitaに投稿
            article_data = {
                "title": title,
                "content": content,
                "tags": tags
            }
            
            result = self.publisher.publish_article(
                article_data,
                private=False,
                tweet=False
            )
            
            if result.get("success"):
                # 状態を更新
                state = self._load_state()
                draft_filename = Path(draft_file).name
                
                if draft_filename not in state.get("published_drafts", []):
                    state.setdefault("published_drafts", []).append(draft_filename)
                
                state["last_published_date"] = datetime.now().isoformat()
                state["current_index"] = state.get("current_index", 0) + 1
                self._save_state(state)
                
                print(f"[SCHEDULER] 投稿成功: {result.get('url')}")
                print(f"[SCHEDULER] 次回の投稿予定: {self.get_next_draft()}")
                
                return {
                    "success": True,
                    "draft_file": draft_file,
                    "title": title,
                    "url": result.get("url"),
                    "item_id": result.get("item_id"),
                    "published": True
                }
            else:
                print(f"[ERROR] 投稿失敗: {result.get('error')}")
                return {
                    "success": False,
                    "draft_file": draft_file,
                    "error": result.get("error")
                }
                
        except Exception as e:
            print(f"[ERROR] ドラフト投稿エラー: {str(e)}")
            import traceback
            traceback.print_exc()
            return {
                "success": False,
                "error": str(e)
            }
    
    def get_status(self) -> Dict[str, Any]:
        """現在の投稿状態を取得"""
        state = self._load_state()
        draft_files = self.get_draft_files()
        
        return {
            "total_drafts": len(draft_files),
            "published_count": len(state.get("published_drafts", [])),
            "current_index": state.get("current_index", 0),
            "next_draft": self.get_next_draft(),
            "last_published_date": state.get("last_published_date"),
            "published_drafts": state.get("published_drafts", [])
        }
