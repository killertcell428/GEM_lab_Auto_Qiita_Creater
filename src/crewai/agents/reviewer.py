"""Reviewer Agent - 品質チェック・校正エージェント"""
from crewai import Agent
from src.crewai.tools.gemini_tool import call_gemini_api
from src.crewai.tools.state_tool import load_article_state
from src.config_loader import get_config


def create_reviewer_agent() -> Agent:
    """
    Reviewer Agentを作成
    
    Returns:
        Reviewer Agent
    """
    config = get_config()
    crewai_config = config.get("crewai", {})
    agent_config = crewai_config.get("agents", {}).get("reviewer", {})
    
    return Agent(
        role="品質チェック・校正者",
        goal="記事の品質をチェックし、改善提案を行う。技術的正確性、構成の論理性、読みやすさを評価する。",
        backstory="""あなたはGEM Lab（遺伝生態モンスター研究所）の品質チェック・校正者です。
        記事の品質を保証し、読者にとって価値のある記事になるよう改善提案を行います。
        
        あなたの役割:
        - 記事の技術的正確性をチェックする
        - 構成の論理性をチェックする
        - 読みやすさ・分かりやすさを評価する
        - 改善提案を生成する
        - コード例の正確性を確認する
        
        チェック項目:
        1. 技術的な内容の正確性
        2. 構成の論理性（導入→説明→実装→まとめの流れ）
        3. コード例の完全性（コピペで動くレベルか）
        4. 文字数（最低3000文字以上）
        5. セクションの完全性（すべてのセクションが含まれているか）
        6. 過去記事との連続性
        7. 伸びている記事の特徴を反映しているか""",
        tools=[
            call_gemini_api,
            load_article_state
        ],
        verbose=agent_config.get("verbose", True),
        allow_delegation=False,
        max_iter=agent_config.get("max_iter", 3)
    )

