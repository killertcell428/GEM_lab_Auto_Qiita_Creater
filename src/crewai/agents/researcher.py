"""Researcher Agent - 技術調査員エージェント"""
from crewai import Agent
from src.crewai.tools.gemini_tool import call_gemini_api
from src.crewai.tools.qiita_tool import fetch_qiita_articles, fetch_trending_articles
from src.crewai.tools.state_tool import load_article_state
from src.config_loader import get_config


def create_researcher_agent() -> Agent:
    """
    Researcher Agentを作成
    
    Returns:
        Researcher Agent
    """
    config = get_config()
    crewai_config = config.get("crewai", {})
    agent_config = crewai_config.get("agents", {}).get("researcher", {})
    
    return Agent(
        role="技術調査員",
        goal="指定されたトピックについて深く調査し、記事作成に必要な情報を収集する。最新の技術トレンド、ツール、研究手法を調査し、バイオインフォマティクス領域への応用可能性を探求する。",
        backstory="""あなたはGEM Lab（遺伝生態モンスター研究所）の技術調査員です。
        新大陸に生息する未知生命"モンスター生命体"の遺伝子・生態・進化を解明するため、
        最新の技術トレンドや研究手法を調査します。
        
        あなたの役割:
        - 指定されたトピックについて深く調査する
        - 関連技術・ツールの情報を収集する
        - 最新トレンドを把握する
        - 過去記事との関連性を調査する
        - バイオインフォマティクス領域への応用可能性を探求する
        
        調査結果は、後続のPlanner Agentが記事構成を設計する際に使用されます。
        そのため、実装可能な具体的な情報を含め、技術的な正確性を保つことが重要です。""",
        tools=[
            call_gemini_api,
            fetch_qiita_articles,
            fetch_trending_articles,
            load_article_state
        ],
        verbose=agent_config.get("verbose", True),
        allow_delegation=False,
        max_iter=agent_config.get("max_iter", 3)
    )

