"""Writer Agent - 記事執筆エージェント"""
from crewai import Agent
from src.crewai.tools.gemini_tool import call_gemini_api
from src.crewai.tools.qiita_tool import fetch_qiita_articles
from src.crewai.tools.state_tool import load_article_state
from src.config_loader import get_config


def create_writer_agent() -> Agent:
    """
    Writer Agentを作成
    
    Returns:
        Writer Agent
    """
    config = get_config()
    crewai_config = config.get("crewai", {})
    agent_config = crewai_config.get("agents", {}).get("writer", {})
    
    return Agent(
        role="記事執筆者",
        goal="構成プランに基づいて、Qiita読者にとって価値のある技術記事を執筆する。実践的で再現可能な内容を含める。",
        backstory="""あなたはGEM Lab（遺伝生態モンスター研究所）の記事執筆者です。
        新大陸に生息する未知生命"モンスター生命体"の遺伝子・生態・進化を解析する調査活動を行っています。
        
        あなたの役割:
        - 構成プランに基づいて記事本文を執筆する
        - Markdown形式で出力する
        - コード例・図表の挿入指示を含める
        - GEM Lab世界観を自然に織り込む（技術解説が主軸）
        - 実践的で再現可能な内容にする
        
        重要: この記事はQiitaに投稿する技術記事です。
        技術的な内容が主軸であり、世界観の説明にならないよう注意してください。
        技術にフォーカスし、世界観は自然に織り込む程度にしてください。""",
        tools=[
            call_gemini_api,
            fetch_qiita_articles,
            load_article_state
        ],
        verbose=agent_config.get("verbose", True),
        allow_delegation=False,
        max_iter=agent_config.get("max_iter", 5)
    )

