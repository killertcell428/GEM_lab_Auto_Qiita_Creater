"""Knowledge Extractor Agent - 実装ノウハウ抽出エージェント"""
from crewai import Agent
from src.crewai.tools.gemini_tool import call_gemini_api
from src.crewai.tools.code_execution_tool import execute_python_code
from src.config_loader import get_config


def create_knowledge_extractor_agent() -> Agent:
    """
    Knowledge Extractor Agentを作成
    
    Returns:
        Knowledge Extractor Agent
    """
    config = get_config()
    crewai_config = config.get("crewai", {})
    agent_config = crewai_config.get("agents", {}).get("knowledge_extractor", {})
    
    return Agent(
        role="ノウハウ抽出者",
        goal="実装過程で得た知見やノウハウを抽出し、読者に価値のある情報を提供する。",
        backstory="""あなたはGEM Lab（遺伝生態モンスター研究所）のノウハウ抽出者です。
        実装検証の過程で得られた知見を整理し、読者が実際に使えるノウハウを抽出します。
        
        あなたの役割:
        - 実装過程でつまずいたポイントを特定
        - 最適化のヒントやベストプラクティスを抽出
        - よくあるエラーと対処法を整理
        - 公式ドキュメントに書かれていない実践的な知見を抽出
        - 読者が実際に使える「生きた知識」を提供
        
        抽出したノウハウは、記事の実践編やFAQセクションに反映されます。
        理論的な説明だけでなく、実際に実装して得た知見を重視してください。""",
        tools=[
            call_gemini_api,
            execute_python_code
        ],
        verbose=agent_config.get("verbose", True),
        allow_delegation=False,
        max_iter=agent_config.get("max_iter", 3)
    )

