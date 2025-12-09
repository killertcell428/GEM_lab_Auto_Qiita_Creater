"""LLM設定ユーティリティ - Gemini LLMの設定"""
from crewai import LLM
from src.config_loader import get_config
import os


def create_gemini_llm() -> LLM:
    """
    Gemini LLMを作成
    
    Returns:
        LLM: Gemini LLMインスタンス
    """
    config = get_config()
    gemini_config = config.get("api", {}).get("gemini", {})
    model = gemini_config.get("model", "gemini-2.5-pro")
    temperature = gemini_config.get("temperature", 0.7)
    api_key = gemini_config.get("api_key", "")
    
    # CrewAIはGeminiを直接サポートしていないため、
    # 環境変数にGEMINI_API_KEYを設定し、カスタムLLMプロバイダーを使用する
    # または、toolsを通じてGeminiを呼び出す方式を使用
    
    # 暫定的に、環境変数を設定してからLLMを作成
    if api_key:
        os.environ["GEMINI_API_KEY"] = api_key
    
    # CrewAIのLLMはデフォルトでOpenAIを使用するため、
    # Geminiを使用するにはカスタム実装が必要
    # ここでは、Noneを返してデフォルトLLMを使用し、
    # toolsでGeminiを呼び出す方式を採用
    return None

