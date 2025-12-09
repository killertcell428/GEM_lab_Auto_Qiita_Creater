"""PDCA State管理 - PDCAサイクルの状態を管理するクラス"""
from dataclasses import dataclass, asdict
from typing import Optional, Dict, Any, List
from datetime import datetime
import json
from pathlib import Path
from src.config_loader import get_config


@dataclass
class PDCAState:
    """PDCAサイクルの状態を管理するクラス"""
    cycle_id: str
    article_id: str
    plan_result: Optional[Dict[str, Any]] = None
    do_result: Optional[Dict[str, Any]] = None
    check_result: Optional[Dict[str, Any]] = None
    act_result: Optional[Dict[str, Any]] = None
    improvements: Optional[List[str]] = None
    created_at: Optional[str] = None
    updated_at: Optional[str] = None
    completed_at: Optional[str] = None
    
    def __post_init__(self):
        if self.created_at is None:
            self.created_at = datetime.now().isoformat()
        if self.updated_at is None:
            self.updated_at = datetime.now().isoformat()
        if self.improvements is None:
            self.improvements = []
    
    def to_dict(self) -> Dict[str, Any]:
        """Stateを辞書形式に変換"""
        return asdict(self)
    
    def save(self, base_dir: Optional[Path] = None):
        """StateをJSONファイルに保存"""
        if base_dir is None:
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            pdca_dir = state_config.get("pdca_dir", "data/state/pdca")
            base_dir = Path(pdca_dir)
        
        base_dir.mkdir(parents=True, exist_ok=True)
        file_path = base_dir / f"{self.cycle_id}.json"
        self.updated_at = datetime.now().isoformat()
        
        with open(file_path, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, ensure_ascii=False, indent=2)
    
    @classmethod
    def load(cls, cycle_id: str, base_dir: Optional[Path] = None) -> Optional["PDCAState"]:
        """StateをJSONファイルから読み込み"""
        if base_dir is None:
            config = get_config()
            crewai_config = config.get("crewai", {})
            state_config = crewai_config.get("state", {})
            pdca_dir = state_config.get("pdca_dir", "data/state/pdca")
            base_dir = Path(pdca_dir)
        
        file_path = base_dir / f"{cycle_id}.json"
        if not file_path.exists():
            return None
        
        with open(file_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        return cls(**data)
    
    def update_plan(self, plan_result: Dict[str, Any]):
        """Plan結果を更新"""
        self.plan_result = plan_result
        self.updated_at = datetime.now().isoformat()
    
    def update_do(self, do_result: Dict[str, Any]):
        """Do結果を更新"""
        self.do_result = do_result
        self.updated_at = datetime.now().isoformat()
    
    def update_check(self, check_result: Dict[str, Any]):
        """Check結果を更新"""
        self.check_result = check_result
        self.updated_at = datetime.now().isoformat()
    
    def update_act(self, act_result: Dict[str, Any]):
        """Act結果を更新"""
        self.act_result = act_result
        self.updated_at = datetime.now().isoformat()
        self.completed_at = datetime.now().isoformat()
    
    def add_improvement(self, improvement: str):
        """改善点を追加"""
        if self.improvements is None:
            self.improvements = []
        self.improvements.append(improvement)
        self.updated_at = datetime.now().isoformat()

