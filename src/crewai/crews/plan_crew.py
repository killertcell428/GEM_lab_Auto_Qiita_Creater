"""Plan Crew - Plan PhaseのCrew実装"""
from crewai import Crew, Process
from typing import Dict, Any, Optional
from src.crewai.tasks.plan_tasks import create_research_task, create_planning_task
from src.crewai.tasks.verification_tasks import create_web_research_task
from src.crewai.state.article_state import ArticleState
from src.config_loader import get_config
import json


def execute_plan_phase(article_state: ArticleState, topic: str, context: Optional[Dict[str, Any]] = None, draft_file: Optional[str] = None) -> ArticleState:
    """
    Plan Phaseを実行
    
    Args:
        article_state: Article State
        topic: 調査するトピック
        context: 追加のコンテキスト情報
    
    Returns:
        更新されたArticle State
    """
    if context is None:
        context = {}
    
    # ドラフトファイルを読み込む（指定されている場合）
    if draft_file:
        try:
            from src.generator.draft_loader import load_draft_from_file
            draft_data = load_draft_from_file(draft_file)
            context["draft_content"] = draft_data.get("raw_content", draft_data.get("content", ""))
            context["draft_title"] = draft_data.get("title")
            context["draft_tags"] = draft_data.get("tags", [])
            print(f"[PLAN] ドラフトファイルを読み込みました: {draft_file}")
        except Exception as e:
            print(f"[WARN] ドラフトファイルの読み込みに失敗: {str(e)}")
            context["draft_content"] = ""
    else:
        context["draft_content"] = ""
    
    # 過去記事データを取得
    try:
        from src.storage.qiita_items_manager import QiitaItemsManager
        qiita_manager = QiitaItemsManager()
        past_articles_summary = qiita_manager.get_items_summary(limit=10)
        context["past_articles"] = past_articles_summary
    except Exception as e:
        print(f"[WARN] 過去記事の取得をスキップ: {str(e)}")
        context["past_articles"] = "なし"
    
    # 伸びている記事の特徴を取得
    try:
        from src.analyzer.article_analyzer import ArticleAnalyzer
        analyzer = ArticleAnalyzer()
        features = analyzer.load_features()
        if features:
            context["trending_features"] = json.dumps(features, ensure_ascii=False)
        else:
            context["trending_features"] = "（分析データなし）"
    except Exception as e:
        print(f"[WARN] トレンド特徴の取得をスキップ: {str(e)}")
        context["trending_features"] = "（分析データなし）"
    
    # 記事タイプを設定
    context.setdefault("article_type", "一般技術記事")
    context.setdefault("article_type_description", "")
    
    # タスク作成
    # 1. Web検索タスク（公式ドキュメント・実装例の検索）
    draft_content = context.get("draft_content", "")
    web_research_task = create_web_research_task(topic, draft_content, context)
    
    # 2. 技術調査タスク（従来の調査）
    research_task = create_research_task(topic, context)
    
    # 3. 構成設計タスク（Web検索結果と調査結果を使用）
    planning_task = create_planning_task("", context)
    
    # Crew作成
    config = get_config()
    crewai_config = config.get("crewai", {})
    crew_config = crewai_config.get("crews", {}).get("plan_crew", {})
    
    crew = Crew(
        agents=[web_research_task.agent, research_task.agent, planning_task.agent],
        tasks=[web_research_task, research_task, planning_task],
        process=Process.sequential,
        verbose=crew_config.get("verbose", True)
    )
    
    # 実行
    print("[PLAN] Plan Phaseを開始...")
    result = crew.kickoff()
    
    # 結果をパース
    # CrewAIの結果は通常、最後のタスクの出力を含む
    # 各タスクの出力を個別に取得
    web_research_result = ""
    research_report = ""
    planning_result = str(result)  # 最後のタスクの出力
    
    # Web検索結果を取得
    if len(crew.tasks) > 0 and hasattr(crew.tasks[0], 'output') and crew.tasks[0].output:
        web_research_result = str(crew.tasks[0].output.raw) if hasattr(crew.tasks[0].output, 'raw') else str(crew.tasks[0].output)
    
    # 技術調査結果を取得
    if len(crew.tasks) > 1 and hasattr(crew.tasks[1], 'output') and crew.tasks[1].output:
        research_report = str(crew.tasks[1].output.raw) if hasattr(crew.tasks[1].output, 'raw') else str(crew.tasks[1].output)
    
    # Web検索結果をcontextに保存（Do Phaseで使用）
    context["web_research_results"] = web_research_result
    
    # JSONを抽出（planning_resultから）
    try:
        # JSON部分を抽出
        if "{" in planning_result and "}" in planning_result:
            json_start = planning_result.find("{")
            json_end = planning_result.rfind("}") + 1
            plan_json = planning_result[json_start:json_end]
            plan = json.loads(plan_json)
        else:
            plan = {"title": topic, "sections": []}
    except Exception as e:
        print(f"[WARN] 計画結果のパースに失敗: {str(e)}")
        plan = {"title": topic, "sections": []}
    
    # State更新
    article_state.topic = topic
    article_state.research_report = research_report
    article_state.plan = plan
    article_state.update_phase("do")
    article_state.save()
    
    print(f"[PLAN] Plan Phase完了: {article_state.article_id}")
    
    return article_state

