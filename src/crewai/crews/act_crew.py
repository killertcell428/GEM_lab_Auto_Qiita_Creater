"""Act Crew - Act PhaseのCrew実装"""
from crewai import Crew, Process
from typing import Dict, Any, Optional
from src.crewai.tasks.act_tasks import create_improvement_abstraction_task, create_plan_reflection_task
from src.crewai.state.article_state import ArticleState
from src.crewai.state.pdca_state import PDCAState
from src.config_loader import get_config
import json
from datetime import datetime


def execute_act_phase(
    article_state: ArticleState,
    pdca_state: PDCAState,
    context: Optional[Dict[str, Any]] = None
) -> tuple:
    """
    Act Phaseを実行
    
    Args:
        article_state: Article State（Check Phase完了後の状態）
        pdca_state: PDCA State
        context: 追加のコンテキスト情報
    
    Returns:
        更新されたArticle StateとPDCA State
    """
    if context is None:
        context = {}
    
    # 次回トピック候補を取得（Domain Trend Analysis結果から）
    analysis_results = article_state.analysis_results or {}
    trend_analysis = analysis_results.get("trend", {})
    next_topic_candidates = trend_analysis.get("next_topic_candidates", [])
    context["next_topic_candidates"] = next_topic_candidates
    
    # タスク作成
    improvement_task = create_improvement_abstraction_task(article_state, pdca_state, context)
    reflection_task = create_plan_reflection_task({}, context)  # 後でimprovement_taskの結果を使用
    
    # Crew作成
    config = get_config()
    crewai_config = config.get("crewai", {})
    crew_config = crewai_config.get("crews", {}).get("act_crew", {})
    
    crew = Crew(
        agents=[improvement_task.agent, reflection_task.agent],
        tasks=[improvement_task, reflection_task],
        process=Process.sequential,
        verbose=crew_config.get("verbose", True)
    )
    
    # 実行
    print("[ACT] Act Phaseを開始...")
    result = crew.kickoff()
    
    # 結果をパース
    result_str = str(result)
    improvement_guidelines = {}
    reflection_result = {}
    
    try:
        # 改善指針を抽出
        if "{" in result_str and "}" in result_str:
            json_start = result_str.find("{")
            json_end = result_str.rfind("}") + 1
            improvement_guidelines = json.loads(result_str[json_start:json_end])
    except Exception as e:
        print(f"[WARN] 改善指針のパースに失敗: {str(e)}")
    
    # Lessons Learnedを抽出
    lessons_learned = improvement_guidelines.get("lessons_learned", [])
    for lesson in lessons_learned:
        article_state.add_lesson_learned(lesson)
    
    # State更新
    article_state.update_phase("completed")
    article_state.save()
    
    # PDCA State更新
    pdca_state.update_act({
        "improvement_guidelines": improvement_guidelines,
        "reflection_result": reflection_result,
        "lessons_learned": lessons_learned
    })
    pdca_state.save()
    
    print(f"[ACT] Act Phase完了: {article_state.article_id}")
    
    return article_state, pdca_state

