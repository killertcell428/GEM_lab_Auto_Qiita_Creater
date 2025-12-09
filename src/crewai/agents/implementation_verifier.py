"""Implementation Verifier Agent - コード実装検証エージェント"""
from crewai import Agent
from src.crewai.tools.code_execution_tool import execute_python_code, verify_code_implementation
from src.crewai.tools.gemini_tool import call_gemini_api
from src.config_loader import get_config


def create_implementation_verifier_agent() -> Agent:
    """
    Implementation Verifier Agentを作成
    
    Returns:
        Implementation Verifier Agent
    """
    config = get_config()
    crewai_config = config.get("crewai", {})
    agent_config = crewai_config.get("agents", {}).get("implementation_verifier", {})
    
    return Agent(
        role="実装検証者",
        goal="実際にコードを実行して動作確認し、実装ノウハウを抽出する。",
        backstory="""あなたはGEM Lab（遺伝生態モンスター研究所）の実装検証者です。
        記事に含めるコードが実際に動作することを確認し、実装過程で得た知見を抽出します。
        
        あなたの役割:
        - コードを実際に実行して動作確認
        - エラーが発生した場合、原因を特定して修正案を提示
        - 実行時間やリソース使用量を確認
        - 実装過程で得たノウハウ（つまずきポイント、最適化のヒントなど）を抽出
        - よくあるエラーと対処法を記録
        
        検証結果は、記事の実践編に反映され、読者が実際に使えるコードを保証します。
        エラーが発生した場合は、必ず修正案を提示してください。""",
        tools=[
            execute_python_code,
            verify_code_implementation,
            call_gemini_api
        ],
        verbose=agent_config.get("verbose", True),
        allow_delegation=False,
        max_iter=agent_config.get("max_iter", 5)
    )

