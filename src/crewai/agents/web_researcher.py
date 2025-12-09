"""Web Researcher Agent - 公式ドキュメント・ブログ検索エージェント"""
from crewai import Agent
from src.crewai.tools.web_search_tool import search_web, search_official_documentation, search_implementation_examples
from src.crewai.tools.gemini_tool import call_gemini_api
from src.config_loader import get_config


def create_web_researcher_agent() -> Agent:
    """
    Web Researcher Agentを作成
    
    Returns:
        Web Researcher Agent
    """
    config = get_config()
    crewai_config = config.get("crewai", {})
    agent_config = crewai_config.get("agents", {}).get("web_researcher", {})
    
    return Agent(
        role="Web調査員",
        goal="公式ドキュメント、ブログ記事、GitHubリポジトリなどを検索し、最新の技術情報と実装例を収集する。",
        backstory="""あなたはGEM Lab（遺伝生態モンスター研究所）のWeb調査員です。
        技術記事の品質を向上させるため、公式ドキュメントや実装例を徹底的に調査します。
        
        あなたの役割:
        - 公式ドキュメントを検索して最新の情報を取得
        - 実装例やブログ記事を検索してベストプラクティスを収集
        - GitHubリポジトリを検索して実際のコード例を探す
        - 技術的な正確性を保つための情報源を特定
        
        調査結果は、後続の実装検証と記事執筆に使用されます。
        そのため、信頼性の高い情報源（公式ドキュメント、有名ブログ、GitHubのスター数が多いリポジトリなど）を優先してください。""",
        tools=[
            search_web,
            search_official_documentation,
            search_implementation_examples,
            call_gemini_api
        ],
        verbose=agent_config.get("verbose", True),
        allow_delegation=False,
        max_iter=agent_config.get("max_iter", 5)
    )

