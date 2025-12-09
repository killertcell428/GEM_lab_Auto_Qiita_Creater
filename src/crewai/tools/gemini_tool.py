"""Gemini API呼び出しツール - CrewAI Toolとして実装"""
try:
    from crewai_tools import tool
except ImportError:
    from crewai.tools import tool
from typing import Dict, Any, Optional
import google.generativeai as genai
from src.config_loader import get_config


@tool
def call_gemini_api(
    prompt: str,
    model: Optional[str] = None,
    max_tokens: Optional[int] = None,
    temperature: Optional[float] = None
) -> str:
    """
    Gemini APIを呼び出してテキストを生成する
    
    Args:
        prompt: プロンプト
        model: 使用するモデル名（Noneの場合は設定ファイルから取得）
        max_tokens: 最大トークン数（Noneの場合は設定ファイルから取得）
        temperature: 温度パラメータ（Noneの場合は設定ファイルから取得）
    
    Returns:
        生成されたテキスト
    """
    config = get_config()
    gemini_config = config.get("api", {}).get("gemini", {})
    api_key = gemini_config.get("api_key", "")
    
    if not api_key:
        raise ValueError("GEMINI_API_KEYが設定されていません")
    
    genai.configure(api_key=api_key)
    
    # デフォルト値の取得
    if model is None:
        model = gemini_config.get("model", "gemini-2.5-pro")
    if max_tokens is None:
        max_tokens = gemini_config.get("max_tokens", 200000)
    if temperature is None:
        temperature = gemini_config.get("temperature", 0.7)
    
    # モデル名の正規化
    if model.startswith("models/"):
        model = model.replace("models/", "")
    
    try:
        model_instance = genai.GenerativeModel(model)
        generation_config = genai.types.GenerationConfig(
            temperature=temperature,
            max_output_tokens=max_tokens,
        )
        
        response = model_instance.generate_content(prompt, generation_config=generation_config)
        return response.text
    except Exception as e:
        raise Exception(f"Gemini API呼び出しに失敗しました: {str(e)}")

