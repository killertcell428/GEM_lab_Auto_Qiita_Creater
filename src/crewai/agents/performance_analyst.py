"""Performance Analyst Agent - パフォーマンス分析エージェント"""
from crewai import Agent
from src.crewai.tools.gemini_tool import call_gemini_api
from src.crewai.tools.qiita_tool import fetch_qiita_articles
from src.crewai.tools.state_tool import load_article_state
from src.config_loader import get_config


def create_performance_analyst_agent() -> Agent:
    """
    Performance Analyst Agentを作成
    
    Returns:
        Performance Analyst Agent
    """
    config = get_config()
    crewai_config = config.get("crewai", {})
    agent_config = crewai_config.get("agents", {}).get("performance_analyst", {})
    
    return Agent(
        role="パフォーマンス分析者",
        goal="投稿後の記事のKPIを分析し、高パフォーマンス記事の特徴を抽出する。",
        backstory="""あなたはGEM Lab（遺伝生態モンスター研究所）のパフォーマンス分析者です。
        投稿後の記事のパフォーマンス（いいね数、閲覧数、コメント数）を分析し、
        高パフォーマンス記事の特徴を抽出して、今後の記事改善に役立てます。
        
        あなたの役割:
        - いいね数・閲覧数・コメント数を分析する
        - 時系列でのパフォーマンス追跡
        - 高パフォーマンス記事の特徴抽出
        - 低パフォーマンス記事の原因分析
        - 改善のヒントを提案する""",
        tools=[
            call_gemini_api,
            fetch_qiita_articles,
            load_article_state
        ],
        verbose=agent_config.get("verbose", True),
        allow_delegation=False,
        max_iter=agent_config.get("max_iter", 3)
    )

