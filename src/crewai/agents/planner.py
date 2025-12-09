"""Planner Agent - 記事構成設計エージェント"""
from crewai import Agent
from src.crewai.tools.gemini_tool import call_gemini_api
from src.crewai.tools.qiita_tool import fetch_qiita_articles
from src.crewai.tools.state_tool import load_article_state
from src.config_loader import get_config


def create_planner_agent() -> Agent:
    """
    Planner Agentを作成
    
    Returns:
        Planner Agent
    """
    config = get_config()
    crewai_config = config.get("crewai", {})
    agent_config = crewai_config.get("agents", {}).get("planner", {})
    
    return Agent(
        role="記事構成設計者",
        goal="調査結果を基に、Qiita記事の構成を設計する。読者が理解しやすく、実践的な記事構成を作成する。",
        backstory="""あなたはGEM Lab（遺伝生態モンスター研究所）の記事構成設計者です。
        技術調査員が収集した情報を基に、Qiita読者にとって価値のある記事構成を設計します。
        
        あなたの役割:
        - 調査結果を基に記事構成を設計する
        - セクション分割・見出しを設計する
        - 読者ターゲットを明確化する
        - 過去記事との連続性を設計する
        - コード例の配置を計画する
        
        設計する記事構成は、後続のWriter Agentが記事を執筆する際の指針となります。
        そのため、具体的で実装可能な構成を提案することが重要です。""",
        tools=[
            call_gemini_api,
            fetch_qiita_articles,
            load_article_state
        ],
        verbose=agent_config.get("verbose", True),
        allow_delegation=False,
        max_iter=agent_config.get("max_iter", 3)
    )

