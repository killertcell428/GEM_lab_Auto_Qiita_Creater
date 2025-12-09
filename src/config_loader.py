"""設定ファイルと環境変数の読み込みユーティリティ"""
import json
import os
from pathlib import Path
from typing import Dict, Any
from dotenv import load_dotenv


def load_config() -> Dict[str, Any]:
    """設定ファイルと環境変数を読み込む"""
    # 環境変数を読み込む
    env_path = Path(".env")
    if env_path.exists():
        load_dotenv(env_path)
    
    # 設定ファイルを読み込む
    config_path = Path("config/config.json")
    if not config_path.exists():
        # 例示ファイルからコピーするよう促す
        example_path = Path("config/config.example.json")
        if example_path.exists():
            raise FileNotFoundError(
                f"設定ファイルが見つかりません。{example_path} を config/config.json にコピーして設定してください。"
            )
        else:
            raise FileNotFoundError("設定ファイルが見つかりません。")
    
    with open(config_path, "r", encoding="utf-8") as f:
        config = json.load(f)
    
    # 環境変数からAPIキーを取得
    config["api"]["gemini"]["api_key"] = os.getenv("GEMINI_API_KEY", "")
    config["api"]["qiita"]["access_token"] = os.getenv("QIITA_ACCESS_TOKEN", "")
    
    # カスタムリサーチAPI設定
    custom_key = os.getenv("CUSTOM_RESEARCH_API_KEY", "")
    custom_url = os.getenv("CUSTOM_RESEARCH_API_URL", "")
    if custom_key:
        config["api"]["custom_research"]["api_key"] = custom_key
    if custom_url:
        config["api"]["custom_research"]["endpoint"] = custom_url
    
    return config


def get_config() -> Dict[str, Any]:
    """設定を取得（シングルトン的な使い方）"""
    if not hasattr(get_config, "_config"):
        get_config._config = load_config()
    return get_config._config

