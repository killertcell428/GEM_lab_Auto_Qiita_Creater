"""Domain Trend Analyst Agent - ドメイントレンド分析エージェント"""
from crewai import Agent
from src.crewai.tools.gemini_tool import call_gemini_api
from src.crewai.tools.qiita_tool import fetch_trending_articles, fetch_qiita_articles
from src.crewai.tools.state_tool import load_article_state
from src.config_loader import get_config


def create_domain_trend_analyst_agent() -> Agent:
    """
    Domain Trend Analyst Agentを作成
    
    Returns:
        Domain Trend Analyst Agent
    """
    config = get_config()
    crewai_config = config.get("crewai", {})
    agent_config = crewai_config.get("agents", {}).get("domain_trend_analyst", {})
    
    return Agent(
        role="ドメイントレンド分析者",
        goal="バイオインフォマティクス領域のトレンドを把握し、次回記事のトピック候補を提案する。",
        backstory="""あなたはGEM Lab（遺伝生態モンスター研究所）のドメイントレンド分析者です。
        バイオインフォマティクス領域の技術トレンドを分析し、
        次回記事のトピック候補を提案します。
        
        あなたの役割:
        - バイオインフォマティクス領域のトレンド把握
        - 技術トレンドの時系列分析
        - 次回記事のトピック候補提案
        - ドメイン知識の更新
        - 読者のニーズの把握""",
        tools=[
            call_gemini_api,
            fetch_trending_articles,
            fetch_qiita_articles,
            load_article_state
        ],
        verbose=agent_config.get("verbose", True),
        allow_delegation=False,
        max_iter=agent_config.get("max_iter", 3)
    )

