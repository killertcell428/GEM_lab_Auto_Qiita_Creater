"""Human Feedback Queue管理 - 人間のフィードバックを管理するクラス"""
from dataclasses import dataclass, asdict
from typing import List, Dict, Any, Optional
from datetime import datetime
import json
from pathlib import Path
from src.config_loader import get_config


@dataclass
class HumanFeedback:
    """人間のフィードバックを管理"""
    feedback_id: str
    article_id: str
    phase: str  # "plan", "do", "check", "act"
    feedback_type: str  # "instruction", "correction", "approval", "rejection"
    content: str  # 自然言語の指示・修正内容
    target_agent: Optional[str] = None  # 対象Agent（Noneの場合はEditor-in-Chiefが判定）
    priority: int = 5  # 1-10（10が最高優先度）
    status: str = "pending"  # "pending", "processing", "completed", "rejected"
    created_at: Optional[str] = None
    processed_at: Optional[str] = None
    task_definition: Optional[Dict[str, Any]] = None  # 変換されたタスク定義
    
    def __post_init__(self):
        if self.created_at is None:
            self.created_at = datetime.now().isoformat()
    
    def to_dict(self) -> Dict[str, Any]:
        """Feedbackを辞書形式に変換"""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "HumanFeedback":
        """辞書からFeedbackを作成"""
        return cls(**data)


class FeedbackQueue:
    """Human Feedbackのキュー管理"""
    
    def __init__(self, base_dir: Optional[Path] = None):
        """初期化"""
        if base_dir is None:
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            feedback_dir = state_config.get("feedback_dir", "data/state/feedback")
            base_dir = Path(feedback_dir)
        
        self.base_dir = base_dir
        self.base_dir.mkdir(parents=True, exist_ok=True)
    
    def add_feedback(self, feedback: HumanFeedback):
        """フィードバックを追加"""
        file_path = self.base_dir / f"{feedback.feedback_id}.json"
        with open(file_path, "w", encoding="utf-8") as f:
            json.dump(feedback.to_dict(), f, ensure_ascii=False, indent=2)
    
    def get_feedback(self, feedback_id: str) -> Optional[HumanFeedback]:
        """フィードバックを取得"""
        file_path = self.base_dir / f"{feedback_id}.json"
        if not file_path.exists():
            return None
        
        with open(file_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        return HumanFeedback.from_dict(data)
    
    def get_pending_feedbacks(self, article_id: Optional[str] = None) -> List[HumanFeedback]:
        """待機中のフィードバックを取得"""
        feedbacks = []
        
        for file_path in self.base_dir.glob("*.json"):
            try:
                with open(file_path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                    if data.get("status") == "pending":
                        if article_id is None or data.get("article_id") == article_id:
                            feedbacks.append(HumanFeedback.from_dict(data))
            except Exception as e:
                print(f"[WARN] フィードバックファイルの読み込みに失敗: {file_path}, {str(e)}")
        
        # 優先度順にソート
        return sorted(feedbacks, key=lambda x: x.priority, reverse=True)
    
    def update_feedback(self, feedback: HumanFeedback):
        """フィードバックを更新"""
        self.add_feedback(feedback)
    
    def mark_as_processing(self, feedback_id: str):
        """フィードバックを処理中にマーク"""
        feedback = self.get_feedback(feedback_id)
        if feedback:
            feedback.status = "processing"
            self.update_feedback(feedback)
    
    def mark_as_completed(self, feedback_id: str, task_definition: Optional[Dict[str, Any]] = None):
        """フィードバックを完了にマーク"""
        feedback = self.get_feedback(feedback_id)
        if feedback:
            feedback.status = "completed"
            feedback.processed_at = datetime.now().isoformat()
            if task_definition:
                feedback.task_definition = task_definition
            self.update_feedback(feedback)
    
    def mark_as_rejected(self, feedback_id: str, reason: Optional[str] = None):
        """フィードバックを却下にマーク"""
        feedback = self.get_feedback(feedback_id)
        if feedback:
            feedback.status = "rejected"
            feedback.processed_at = datetime.now().isoformat()
            if reason:
                feedback.content = f"{feedback.content} [却下理由: {reason}]"
            self.update_feedback(feedback)

