"""Act Phase Tasks - 改善指針の抽象化・反映タスク"""
from crewai import Task
from typing import Dict, Any
from src.crewai.agents.editor_in_chief import create_editor_in_chief_agent
from src.crewai.agents.planner import create_planner_agent
from src.crewai.state.article_state import ArticleState
from src.crewai.state.pdca_state import PDCAState
import json


def create_improvement_abstraction_task(
    article_state: ArticleState,
    pdca_state: "PDCAState",
    context: Dict[str, Any]
) -> Task:
    """
    改善指針の抽象化タスクを作成
    
    Args:
        article_state: Article State
        pdca_state: PDCA State
        context: コンテキスト情報
    
    Returns:
        Improvement Abstraction Task
    """
    editor = create_editor_in_chief_agent()
    
    analysis_results = article_state.analysis_results or {}
    review_result = article_state.review_result or {}
    
    return Task(
        description=f"""
        以下の分析結果を基に、改善指針を抽象化してください：
        
        記事ID: {article_state.article_id}
        
        パフォーマンス分析結果:
        {json.dumps(analysis_results.get('performance', {}), ensure_ascii=False, indent=2)}
        
        外部ベンチマーク分析結果:
        {json.dumps(analysis_results.get('benchmark', {}), ensure_ascii=False, indent=2)}
        
        ドメイントレンド分析結果:
        {json.dumps(analysis_results.get('trend', {}), ensure_ascii=False, indent=2)}
        
        レビュー結果:
        {json.dumps(review_result, ensure_ascii=False, indent=2)}
        
        改善指針の抽象化項目:
        1. 高パフォーマンス記事の共通特徴
        2. 低パフォーマンス記事の共通原因
        3. 次回記事への具体的な改善提案
        4. プロンプト改善の提案
        5. プロセス改善の提案
        
        改善指針はJSON形式で出力してください：
        {{
            "high_performance_patterns": ["パターン1", "パターン2", ...],
            "low_performance_patterns": ["パターン1", "パターン2", ...],
            "next_article_improvements": ["改善1", "改善2", ...],
            "prompt_improvements": ["プロンプト改善1", "改善2", ...],
            "process_improvements": ["プロセス改善1", "改善2", ...],
            "lessons_learned": ["学び1", "学び2", ...]
        }}
        """,
        agent=editor,
        expected_output="改善指針（JSON形式）。高パフォーマンスパターン、低パフォーマンスパターン、次回記事改善、プロンプト改善、プロセス改善、Lessons Learnedを含む。"
    )


def create_plan_reflection_task(
    improvement_guidelines: Dict[str, Any],
    context: Dict[str, Any]
) -> Task:
    """
    次回計画への反映タスクを作成
    
    Args:
        improvement_guidelines: 改善指針
        context: コンテキスト情報
    
    Returns:
        Plan Reflection Task
    """
    planner = create_planner_agent()
    
    return Task(
        description=f"""
        以下の改善指針を次回記事の計画に反映してください：
        
        改善指針:
        {json.dumps(improvement_guidelines, ensure_ascii=False, indent=2)}
        
        次回トピック候補:
        {json.dumps(context.get('next_topic_candidates', []), ensure_ascii=False, indent=2)}
        
        反映項目:
        1. 改善指針を次回記事の構成設計に反映
        2. 高パフォーマンスパターンを活用
        3. 低パフォーマンスパターンを回避
        4. プロンプト改善を反映
        
        反映結果はJSON形式で出力してください：
        {{
            "reflected_improvements": ["反映した改善1", "改善2", ...],
            "next_plan_guidelines": "次回計画のガイドライン",
            "expected_improvements": ["期待される改善1", "改善2", ...]
        }}
        """,
        agent=planner,
        expected_output="次回計画への反映結果（JSON形式）。反映した改善、次回計画ガイドライン、期待される改善を含む。"
    )

