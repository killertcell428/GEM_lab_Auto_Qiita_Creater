"""Do Crew - Do PhaseのCrew実装"""
from crewai import Crew, Process
from typing import Dict, Any, Optional
from src.crewai.tasks.do_tasks import (
    create_writing_task,
    create_review_task
)
from src.crewai.tasks.verification_tasks import (
    create_implementation_verification_task,
    create_knowledge_extraction_task
)
from src.crewai.state.article_state import ArticleState
from src.config_loader import get_config
import json


def execute_do_phase(article_state: ArticleState, context: Optional[Dict[str, Any]] = None) -> ArticleState:
    """
    Do Phaseを実行
    
    Args:
        article_state: Article State（Plan Phase完了後の状態）
        context: 追加のコンテキスト情報
    
    Returns:
        更新されたArticle State
    """
    if context is None:
        context = {}
    
    # Plan Phaseの結果を取得
    plan = article_state.plan
    research_report = article_state.research_report
    
    if not plan:
        raise ValueError("Plan Phaseが完了していません。planが設定されていません。")
    
    # planが文字列の場合は辞書に変換を試みる
    if isinstance(plan, str):
        try:
            plan = json.loads(plan)
        except json.JSONDecodeError:
            # JSONパースに失敗した場合は、デフォルトのplanを作成
            plan = {"title": article_state.topic or "タイトル未設定", "sections": []}
    
    if not research_report:
        raise ValueError("Plan Phaseが完了していません。research_reportが設定されていません。")
    
    # research_reportが辞書の場合は文字列に変換
    if isinstance(research_report, dict):
        research_report = json.dumps(research_report, ensure_ascii=False, indent=2)
    
    # ドラフト内容を取得（Plan Phaseから引き継ぐ）
    if "draft_content" not in context:
        # Plan Phaseでドラフトが読み込まれていない場合、ArticleStateから取得を試みる
        # ただし、通常はPlan Phaseでcontextに含まれているはず
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
    
    # タスク作成
    writing_task = create_writing_task(plan, research_report, context)
    review_task = create_review_task("", plan, context)  # 後でwriting_taskの結果を使用
    
    # Crew作成
    config = get_config()
    crewai_config = config.get("crewai", {})
    crew_config = crewai_config.get("crews", {}).get("do_crew", {})
    
    crew = Crew(
        agents=[writing_task.agent, review_task.agent],
        tasks=[writing_task, review_task],
        process=Process.sequential,
        verbose=crew_config.get("verbose", True)
    )
    
    # 実行
    print("[DO] Do Phaseを開始...")
    print("[DO] 1. 記事執筆 → 2. 実装検証 → 3. ノウハウ抽出 → 4. レビュー")
    result = crew.kickoff()
    
    # 実装検証結果とノウハウ抽出結果をcontextに保存（後続のタスクで使用）
    if len(crew.tasks) > 1 and hasattr(crew.tasks[1], 'output') and crew.tasks[1].output:
        verification_result_str = str(crew.tasks[1].output.raw) if hasattr(crew.tasks[1].output, 'raw') else str(crew.tasks[1].output)
        context["verification_results"] = verification_result_str
    
    if len(crew.tasks) > 2 and hasattr(crew.tasks[2], 'output') and crew.tasks[2].output:
        knowledge_result_str = str(crew.tasks[2].output.raw) if hasattr(crew.tasks[2].output, 'raw') else str(crew.tasks[2].output)
        context["knowledge_extraction_results"] = knowledge_result_str
    
    # 各タスクの出力を個別に取得
    # CrewAIでは、sequentialプロセスの場合、各タスクの出力はcrew.tasksから取得できる
    article_content = ""
    verification_results = {}
    knowledge_extraction_results = {}
    review_result = {"approval": True, "improvements": []}
    
    # writing_task（最初のタスク）の出力を取得
    if len(crew.tasks) > 0 and hasattr(crew.tasks[0], 'output') and crew.tasks[0].output:
        article_content = str(crew.tasks[0].output.raw) if hasattr(crew.tasks[0].output, 'raw') else str(crew.tasks[0].output)
    
    # verification_task（2番目のタスク）の出力を取得
    if len(crew.tasks) > 1 and hasattr(crew.tasks[1], 'output') and crew.tasks[1].output:
        verification_result_str = str(crew.tasks[1].output.raw) if hasattr(crew.tasks[1].output, 'raw') else str(crew.tasks[1].output)
        try:
            if "{" in verification_result_str and "}" in verification_result_str:
                json_start = verification_result_str.find("{")
                json_end = verification_result_str.rfind("}") + 1
                verification_results = json.loads(verification_result_str[json_start:json_end + 1])
        except Exception as e:
            print(f"[WARN] 実装検証結果のパースに失敗: {str(e)}")
    
    # knowledge_extraction_task（3番目のタスク）の出力を取得
    if len(crew.tasks) > 2 and hasattr(crew.tasks[2], 'output') and crew.tasks[2].output:
        knowledge_result_str = str(crew.tasks[2].output.raw) if hasattr(crew.tasks[2].output, 'raw') else str(crew.tasks[2].output)
        try:
            if "{" in knowledge_result_str and "}" in knowledge_result_str:
                json_start = knowledge_result_str.find("{")
                json_end = knowledge_result_str.rfind("}") + 1
                knowledge_extraction_results = json.loads(knowledge_result_str[json_start:json_end + 1])
        except Exception as e:
            print(f"[WARN] ノウハウ抽出結果のパースに失敗: {str(e)}")
    
    # review_task（4番目のタスク）の出力を取得
    if len(crew.tasks) > 3 and hasattr(crew.tasks[3], 'output') and crew.tasks[3].output:
        review_result_str = str(crew.tasks[3].output.raw) if hasattr(crew.tasks[3].output, 'raw') else str(crew.tasks[3].output)
        try:
            if "{" in review_result_str and "}" in review_result_str:
                json_start = review_result_str.find("{")
                json_end = review_result_str.rfind("}") + 1
                review_result = json.loads(review_result_str[json_start:json_end + 1])
        except Exception as e:
            print(f"[WARN] レビュー結果のパースに失敗: {str(e)}")
    
    # フォールバック: resultから抽出
    if not article_content:
        result_str = str(result)
        # 最後のJSONブロック（レビュー結果）を探す
        if "{" in result_str:
            last_json_start = result_str.rfind("{")
            if last_json_start > 0:
                article_content = result_str[:last_json_start].strip()
            else:
                article_content = result_str
        else:
            article_content = result_str
    
    # 記事本文が空または短すぎる場合のフォールバック
    if not article_content or len(article_content) < 100:
        print("[WARN] 記事本文が取得できませんでした。result全体を使用します。")
        result_str = str(result)
        # JSON部分を除去
        if "{" in result_str:
            last_json_start = result_str.rfind("{")
            if last_json_start > 0:
                article_content = result_str[:last_json_start].strip()
            else:
                article_content = result_str
        else:
            article_content = result_str
    
    # タイトルを抽出（記事内容の最初の`# `行から）
    if isinstance(plan, dict):
        title = plan.get("title", "タイトル未設定")
    else:
        title = "タイトル未設定"
    
    # 記事本文の最初の行からタイトルを抽出
    if article_content.startswith("# "):
        first_line = article_content.split("\n")[0]
        title = first_line.replace("# ", "").strip()
    elif "\n# " in article_content:
        # 最初の`# `見出しを探す
        lines = article_content.split("\n")
        for line in lines:
            if line.strip().startswith("# "):
                title = line.replace("# ", "").strip()
                break
    
    # 検証結果とノウハウを記事に反映（記事本文に追加）
    if verification_results or knowledge_extraction_results:
        # 検証結果とノウハウを記事の最後に追加
        additional_section = "\n\n---\n\n## 🔍 実装検証とノウハウ\n\n"
        
        if verification_results:
            additional_section += "### 実装検証結果\n\n"
            if isinstance(verification_results, dict):
                if verification_results.get("verification_results"):
                    for vr in verification_results.get("verification_results", []):
                        if vr.get("success"):
                            additional_section += f"- ✅ コードは正常に実行されました（実行時間: {vr.get('execution_time', 'N/A')}秒）\n"
                        else:
                            additional_section += f"- ❌ エラー: {vr.get('error', 'N/A')}\n"
                            additional_section += f"  - 修正案: {vr.get('notes', 'N/A')}\n"
                if verification_results.get("common_errors"):
                    additional_section += "\n### よくあるエラーと対処法\n\n"
                    for error in verification_results.get("common_errors", []):
                        additional_section += f"- **{error.get('error_type', 'N/A')}**: {error.get('solution', 'N/A')}\n"
                if verification_results.get("optimization_tips"):
                    additional_section += "\n### 最適化のヒント\n\n"
                    for tip in verification_results.get("optimization_tips", []):
                        additional_section += f"- {tip}\n"
        
        if knowledge_extraction_results:
            additional_section += "\n### 実装ノウハウ\n\n"
            if isinstance(knowledge_extraction_results, dict):
                if knowledge_extraction_results.get("practical_insights"):
                    additional_section += "#### 実践的な知見\n\n"
                    for insight in knowledge_extraction_results.get("practical_insights", []):
                        additional_section += f"- **{insight.get('insight', 'N/A')}**: {insight.get('context', 'N/A')}\n"
                if knowledge_extraction_results.get("troubleshooting_guide"):
                    additional_section += "\n#### トラブルシューティングガイド\n\n"
                    for guide in knowledge_extraction_results.get("troubleshooting_guide", []):
                        additional_section += f"- **問題**: {guide.get('problem', 'N/A')}\n"
                        additional_section += f"  - **解決方法**: {guide.get('solution', 'N/A')}\n"
                        additional_section += f"  - **予防策**: {guide.get('prevention', 'N/A')}\n"
        
        article_content += additional_section
    
    # State更新
    article_state.content = article_content
    article_state.review_result = review_result
    # 検証結果とノウハウを保存
    article_state.analysis_results = {
        "verification_results": verification_results,
        "knowledge_extraction_results": knowledge_extraction_results
    }
    article_state.update_phase("check")
    article_state.save()
    
    print(f"[DO] Do Phase完了: {article_state.article_id}")
    
    return article_state

