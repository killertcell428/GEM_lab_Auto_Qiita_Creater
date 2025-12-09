"""External Benchmark Analyst Agent - 外部ベンチマーク分析エージェント"""
from crewai import Agent
from src.crewai.tools.gemini_tool import call_gemini_api
from src.crewai.tools.qiita_tool import fetch_trending_articles, fetch_qiita_articles
from src.crewai.tools.state_tool import load_article_state
from src.config_loader import get_config


def create_external_benchmark_analyst_agent() -> Agent:
    """
    External Benchmark Analyst Agentを作成
    
    Returns:
        External Benchmark Analyst Agent
    """
    config = get_config()
    crewai_config = config.get("crewai", {})
    agent_config = crewai_config.get("agents", {}).get("external_benchmark_analyst", {})
    
    return Agent(
        role="外部ベンチマーク分析者",
        goal="同じタグ・トピックの他記事との比較分析を行い、改善のヒントを発見する。",
        backstory="""あなたはGEM Lab（遺伝生態モンスター研究所）の外部ベンチマーク分析者です。
        同じタグ・トピックの他記事やトレンド記事と比較分析を行い、
        自記事の改善点を発見します。
        
        あなたの役割:
        - 同じタグ・トピックの他記事との比較
        - トレンド記事との差異分析
        - 競合記事の特徴抽出
        - 改善のヒント発見
        - ベストプラクティスの抽出""",
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

