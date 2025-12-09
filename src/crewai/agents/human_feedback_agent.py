"""Human Feedback Agent - 人間指示解析エージェント"""
from crewai import Agent
from src.crewai.tools.gemini_tool import call_gemini_api
from src.crewai.tools.state_tool import load_article_state
from src.config_loader import get_config
from typing import Dict, Any
import json


def create_human_feedback_agent() -> Agent:
    """
    Human Feedback Agentを作成
    
    Returns:
        Human Feedback Agent
    """
    config = get_config()
    crewai_config = config.get("crewai", {})
    agent_config = crewai_config.get("agents", {}).get("human_feedback_agent", {})
    
    return Agent(
        role="人間指示解析エージェント",
        goal="人間の自然言語指示を解析し、具体的なタスクに変換する。適切なAgentに割り当てるタスクを生成する。",
        backstory="""あなたはGEM Lab（遺伝生態モンスター研究所）の人間指示解析エージェントです。
        人間からの自然言語指示を理解し、適切なAgentに割り当てるタスクに変換します。
        
        あなたの役割:
        - 自然言語指示の解析
        - 指示を具体的なタスクに変換
        - 対象Agentの特定
        - 優先度の判定
        - 必要なコンテキスト情報の抽出
        
        人間の指示は様々な形式で来る可能性があります：
        - 「タイトルを変更して」
        - 「コード例を追加して」
        - 「もっと分かりやすく書いて」
        - 「この部分を削除して」
        
        これらを、具体的なタスク定義に変換してください。""",
        tools=[
            call_gemini_api,
            load_article_state
        ],
        verbose=agent_config.get("verbose", True),
        allow_delegation=False,
        max_iter=agent_config.get("max_iter", 2)
    )


def parse_human_feedback(feedback_content: str, article_state: "ArticleState", phase: str) -> Dict[str, Any]:
    """
    人間のフィードバックをタスクに変換
    
    Args:
        feedback_content: 人間の指示（自然言語）
        article_state: Article State
        phase: 現在のPhase
    
    Returns:
        タスク定義（辞書形式）
    """
    agent = create_human_feedback_agent()
    
    prompt = f"""
    以下の人間の指示を解析し、具体的なタスクに変換してください：
    
    指示: {feedback_content}
    現在のPhase: {phase}
    記事ID: {article_state.article_id}
    記事タイトル: {article_state.plan.get('title', '') if article_state.plan else ''}
    
    以下のJSON形式で出力してください：
    {{
        "target_agent": "Agent名（researcher, planner, writer, reviewer, editor_in_chiefなど）",
        "task_description": "具体的なタスク説明",
        "priority": 優先度（1-10、10が最高優先度）,
        "required_context": "必要なコンテキスト情報",
        "expected_output": "期待される出力"
    }}
    """
    
    try:
        result = call_gemini_api(prompt)
        # JSONをパース
        task_definition = json.loads(result)
        return task_definition
    except Exception as e:
        # パースに失敗した場合のフォールバック
        return {
            "target_agent": "editor_in_chief",
            "task_description": feedback_content,
            "priority": 5,
            "required_context": "記事の現在の状態",
            "expected_output": "指示に従った修正結果"
        }

