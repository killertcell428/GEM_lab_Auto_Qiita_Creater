"""Article State管理 - 記事の状態を管理するクラス"""
from dataclasses import dataclass, asdict
from typing import Optional, Dict, Any, List
from datetime import datetime
import json
from pathlib import Path
from src.config_loader import get_config


@dataclass
class ArticleState:
    """記事の状態を管理するクラス"""
    article_id: str
    topic: Optional[str] = None
    research_report: Optional[str] = None
    plan: Optional[Dict[str, Any]] = None
    content: Optional[str] = None
    review_result: Optional[Dict[str, Any]] = None
    qiita_url: Optional[str] = None
    qiita_item_id: Optional[str] = None
    kpi: Optional[Dict[str, Any]] = None
    analysis_results: Optional[Dict[str, Any]] = None
    lessons_learned: Optional[List[str]] = None
    human_feedback: Optional[List[Dict[str, Any]]] = None
    current_phase: str = "plan"  # "plan", "do", "check", "act", "completed"
    created_at: Optional[str] = None
    updated_at: Optional[str] = None
    # 承認フロー関連フィールド
    pending_approval: bool = False
    approval_deadline: Optional[str] = None  # ISO形式
    approval_status: Optional[str] = None  # "pending", "approved", "auto_published"
    scheduled_publish_date: Optional[str] = None  # ISO形式
    
    def __post_init__(self):
        if self.created_at is None:
            self.created_at = datetime.now().isoformat()
        if self.updated_at is None:
            self.updated_at = datetime.now().isoformat()
        if self.human_feedback is None:
            self.human_feedback = []
        if self.lessons_learned is None:
            self.lessons_learned = []
    
    def to_dict(self) -> Dict[str, Any]:
        """Stateを辞書形式に変換"""
        return asdict(self)
    
    def save(self, base_dir: Optional[Path] = None):
        """StateをJSONファイルに保存"""
        if base_dir is None:
            # configがNoneでも落ちないように防御
            config = get_config() or {}
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            articles_dir = state_config.get("articles_dir", "data/state/articles")
            base_dir = Path(articles_dir)
        
        base_dir.mkdir(parents=True, exist_ok=True)
        file_path = base_dir / f"{self.article_id}.json"
        self.updated_at = datetime.now().isoformat()
        
        with open(file_path, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, ensure_ascii=False, indent=2)
    
    @classmethod
    def load(cls, article_id: str, base_dir: Optional[Path] = None) -> "ArticleState":
        """StateをJSONファイルから読み込み"""
        if base_dir is None:
            # configがNoneでも落ちないように防御
            config = get_config() or {}
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            articles_dir = state_config.get("articles_dir", "data/state/articles")
            base_dir = Path(articles_dir)
        
        file_path = base_dir / f"{article_id}.json"
        if not file_path.exists():
            return cls(article_id=article_id)
        
        with open(file_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        return cls(**data)
    
    def update_phase(self, phase: str):
        """Phaseを更新"""
        self.current_phase = phase
        self.updated_at = datetime.now().isoformat()
    
    def add_human_feedback(self, feedback: Dict[str, Any]):
        """Human Feedbackを追加"""
        if self.human_feedback is None:
            self.human_feedback = []
        self.human_feedback.append({
            **feedback,
            "added_at": datetime.now().isoformat()
        })
        self.updated_at = datetime.now().isoformat()
    
    def add_lesson_learned(self, lesson: str):
        """Lessons Learnedを追加"""
        if self.lessons_learned is None:
            self.lessons_learned = []
        self.lessons_learned.append(lesson)
        self.updated_at = datetime.now().isoformat()

