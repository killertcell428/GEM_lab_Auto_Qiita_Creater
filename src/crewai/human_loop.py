"""Human-in-the-Loop処理 - 人間介入の管理"""
from typing import Dict, Any, Optional
from datetime import datetime
from src.crewai.state.feedback_queue import FeedbackQueue, HumanFeedback
from src.crewai.state.article_state import ArticleState
from src.crewai.agents.human_feedback_agent import parse_human_feedback


class HumanLoop:
    """Human-in-the-Loop処理を管理"""
    
    def __init__(self):
        """初期化"""
        self.feedback_queue = FeedbackQueue()
    
    def add_feedback(
        self,
        article_id: str,
        phase: str,
        content: str,
        feedback_type: str = "instruction",
        target_agent: Optional[str] = None,
        priority: int = 5
    ) -> HumanFeedback:
        """
        フィードバックを追加
        
        Args:
            article_id: 記事ID
            phase: 現在のPhase
            content: フィードバック内容（自然言語）
            feedback_type: フィードバックタイプ（"instruction", "correction", "approval", "rejection"）
            target_agent: 対象Agent（Noneの場合は自動判定）
            priority: 優先度（1-10）
        
        Returns:
            HumanFeedback: 追加されたフィードバック
        """
        feedback = HumanFeedback(
            feedback_id=f"{article_id}_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
            article_id=article_id,
            phase=phase,
            feedback_type=feedback_type,
            content=content,
            target_agent=target_agent,
            priority=priority
        )
        self.feedback_queue.add_feedback(feedback)
        return feedback
    
    def process_feedback(
        self,
        article_state: ArticleState,
        feedback: HumanFeedback
    ) -> tuple[ArticleState, Dict[str, Any]]:
        """
        フィードバックを処理
        
        Args:
            article_state: Article State
            feedback: Human Feedback
        
        Returns:
            tuple: (更新されたArticle State, タスク定義)
        """
        # フィードバックを処理中にマーク
        self.feedback_queue.mark_as_processing(feedback.feedback_id)
        
        # タスクに変換
        task_definition = parse_human_feedback(
            feedback.content,
            article_state,
            feedback.phase
        )
        
        # フィードバックを完了にマーク
        self.feedback_queue.mark_as_completed(feedback.feedback_id, task_definition)
        
        # State更新
        article_state.add_human_feedback({
            "feedback_id": feedback.feedback_id,
            "content": feedback.content,
            "task_definition": task_definition,
            "processed_at": datetime.now().isoformat()
        })
        article_state.save()
        
        return article_state, task_definition
    
    def check_and_process_pending_feedbacks(
        self,
        article_state: ArticleState
    ) -> list[Dict[str, Any]]:
        """
        待機中のフィードバックをチェックして処理
        
        Args:
            article_state: Article State
        
        Returns:
            list: 処理されたタスク定義のリスト
        """
        pending = self.feedback_queue.get_pending_feedbacks(article_state.article_id)
        
        processed_tasks = []
        for feedback in pending:
            article_state, task_definition = self.process_feedback(article_state, feedback)
            processed_tasks.append(task_definition)
        
        return processed_tasks
    
    def get_pending_feedbacks(self, article_id: str) -> list[HumanFeedback]:
        """
        待機中のフィードバックを取得
        
        Args:
            article_id: 記事ID
        
        Returns:
            list: 待機中のフィードバックのリスト
        """
        return self.feedback_queue.get_pending_feedbacks(article_id)

