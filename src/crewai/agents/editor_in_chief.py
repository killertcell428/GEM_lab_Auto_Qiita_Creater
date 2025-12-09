"""Editor-in-Chief Agent - 統合・意思決定エージェント"""
from crewai import Agent
from src.crewai.tools.gemini_tool import call_gemini_api
from src.crewai.tools.state_tool import load_article_state, save_article_state, load_pdca_state, save_pdca_state
from src.crewai.tools.qiita_tool import publish_to_qiita
from src.config_loader import get_config


def create_editor_in_chief_agent() -> Agent:
    """
    Editor-in-Chief Agentを作成
    
    Returns:
        Editor-in-Chief Agent
    """
    config = get_config()
    crewai_config = config.get("crewai", {})
    agent_config = crewai_config.get("agents", {}).get("editor_in_chief", {})
    
    return Agent(
        role="編集長",
        goal="全Agentの統括・調整を行い、PDCAサイクルを管理する。各Phaseの承認・進行管理、Human Feedbackの処理、最終的な意思決定を行う。",
        backstory="""あなたはGEM Lab（遺伝生態モンスター研究所）の編集長です。
        記事生成から分析、改善まで、全プロセスを統括し、PDCAサイクルを回します。
        
        あなたの役割:
        - 全Agentの統括・調整
        - タスクの優先順位決定
        - 各Phaseの承認・進行管理
        - Human Feedbackの処理・割り込みタスクの生成
        - PDCAサイクルの進行管理
        - 最終的な意思決定
        - 改善指針の抽象化
        
        あなたは、各Phaseの結果を評価し、次のアクションを決定します。
        人間からのフィードバックがある場合は、それを優先的に処理します。""",
        tools=[
            call_gemini_api,
            load_article_state,
            save_article_state,
            load_pdca_state,
            save_pdca_state,
            publish_to_qiita
        ],
        verbose=agent_config.get("verbose", True),
        allow_delegation=True,  # Editor-in-Chiefは他のAgentに委任可能
        max_iter=agent_config.get("max_iter", 5)
    )

