"""設定APIルーター"""
from fastapi import APIRouter, HTTPException
from typing import Dict, Any
from src.config_loader import get_config
import json
from pathlib import Path

router = APIRouter()


@router.get("/settings")
async def get_settings():
    """設定取得"""
    try:
        config = get_config()
        return config
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.put("/settings")
async def update_settings(settings: Dict[str, Any]):
    """設定更新"""
    try:
        config_file = Path("config/config.json")
        if not config_file.exists():
            raise HTTPException(status_code=404, detail="設定ファイルが見つかりません")
        
        # 既存の設定を読み込んでマージ
        with open(config_file, "r", encoding="utf-8") as f:
            current_config = json.load(f)
        
        # マージ（深いマージが必要な場合は改善が必要）
        current_config.update(settings)
        
        # 保存
        with open(config_file, "w", encoding="utf-8") as f:
            json.dump(current_config, f, ensure_ascii=False, indent=2)
        
        return {"message": "設定を更新しました", "settings": current_config}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

