"""Check Crew - Check PhaseのCrew実装"""
from crewai import Crew, Process
from typing import Dict, Any, Optional
from src.crewai.tasks.check_tasks import (
    create_performance_analysis_task,
    create_external_benchmark_task,
    create_domain_trend_task
)
from src.crewai.state.article_state import ArticleState
from src.config_loader import get_config
import json


def execute_check_phase(article_state: ArticleState, context: Optional[Dict[str, Any]] = None) -> ArticleState:
    """
    Check Phaseを実行
    
    Args:
        article_state: Article State（Do Phase完了後の状態）
        context: 追加のコンテキスト情報
    
    Returns:
        更新されたArticle State
    """
    if context is None:
        context = {}
    
    # Qiita URLが必要（投稿済みの場合）
    if not article_state.qiita_url:
        print("[WARN] 記事がまだQiitaに投稿されていません。Check Phaseをスキップします。")
        article_state.update_phase("act")
        article_state.save()
        return article_state
    
    # タスク作成（3つのAnalyst Agentを並列実行）
    performance_task = create_performance_analysis_task(article_state, context)
    benchmark_task = create_external_benchmark_task(article_state, context)
    trend_task = create_domain_trend_task(article_state, context)
    
    # Crew作成（sequential processで実行）
    # hierarchicalはmanager_llmが必要なため、sequentialに変更
    config = get_config()
    crewai_config = config.get("crewai", {})
    crew_config = crewai_config.get("crews", {}).get("check_crew", {})
    
    # プロセスタイプを設定
    process_type = crew_config.get("process", "sequential")
    if process_type == "hierarchical":
        process_type = "sequential"  # hierarchicalはmanager_llmが必要なため、sequentialにフォールバック
    if process_type == "parallel":
        process_type = "sequential"  # parallelは存在しないため、sequentialにフォールバック
    
    process_map = {
        "sequential": Process.sequential,
        "hierarchical": Process.hierarchical
    }
    process = process_map.get(process_type, Process.sequential)
    
    crew = Crew(
        agents=[
            performance_task.agent,
            benchmark_task.agent,
            trend_task.agent
        ],
        tasks=[performance_task, benchmark_task, trend_task],
        process=process,
        verbose=crew_config.get("verbose", True)
    )
    
    # 実行
    print("[CHECK] Check Phaseを開始...")
    result = crew.kickoff()
    
    # 結果をパース（各タスクの結果を統合）
    # 実際の実装では、resultから各タスクの出力を適切に抽出する必要がある
    analysis_results = {
        "performance": {},
        "benchmark": {},
        "trend": {}
    }
    
    result_str = str(result)
    
    # 各分析結果を抽出（暫定的な実装）
    try:
        # Performance Analysis結果を抽出
        if "performance" in result_str.lower() or "likes_count" in result_str:
            # JSON部分を抽出
            if "{" in result_str and "}" in result_str:
                json_start = result_str.find("{")
                json_end = result_str.rfind("}") + 1
                analysis_results["performance"] = json.loads(result_str[json_start:json_end])
    except Exception as e:
        print(f"[WARN] パフォーマンス分析結果のパースに失敗: {str(e)}")
    
    # State更新
    article_state.analysis_results = analysis_results
    article_state.update_phase("act")
    article_state.save()
    
    print(f"[CHECK] Check Phase完了: {article_state.article_id}")
    
    return article_state

