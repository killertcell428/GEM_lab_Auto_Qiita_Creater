"""プロンプト管理モジュール - data/prompts/からプロンプトファイルを読み込む"""
from pathlib import Path
from typing import Optional, Dict


class PromptLoader:
    """プロンプトファイルを読み込むクラス"""
    
    def __init__(self, prompts_dir: str = "data/prompts"):
        """
        Args:
            prompts_dir: プロンプトファイルを保存するディレクトリ
        """
        self.prompts_dir = Path(prompts_dir)
        self.prompts_dir.mkdir(parents=True, exist_ok=True)
        self._cache: Dict[str, str] = {}
    
    def load_prompt(self, filename: str) -> str:
        """
        プロンプトファイルを読み込む
        
        Args:
            filename: プロンプトファイル名（例: "research_prompt.md"）
            
        Returns:
            str: プロンプトの内容
            
        Raises:
            FileNotFoundError: ファイルが見つからない場合
        """
        # キャッシュをチェック
        if filename in self._cache:
            return self._cache[filename]
        
        # ファイルパスを構築
        prompt_path = self.prompts_dir / filename
        
        if not prompt_path.exists():
            raise FileNotFoundError(
                f"プロンプトファイルが見つかりません: {prompt_path}\n"
                f"data/prompts/ ディレクトリに {filename} を作成してください。"
            )
        
        # ファイルを読み込む
        with open(prompt_path, "r", encoding="utf-8") as f:
            content = f.read().strip()
        
        # キャッシュに保存
        self._cache[filename] = content
        
        return content
    
    def load_prompt_or_default(self, filename: str, default: str) -> str:
        """
        プロンプトファイルを読み込む（存在しない場合はデフォルト値を返す）
        
        Args:
            filename: プロンプトファイル名
            default: デフォルトのプロンプト内容
            
        Returns:
            str: プロンプトの内容
        """
        try:
            return self.load_prompt(filename)
        except FileNotFoundError:
            return default
    
    def format_prompt(self, filename: str, **kwargs) -> str:
        """
        プロンプトファイルを読み込んで変数を置換
        
        Args:
            filename: プロンプトファイル名
            **kwargs: 置換する変数（例: topic="Python", count=5）
            
        Returns:
            str: フォーマットされたプロンプト
        """
        prompt = self.load_prompt(filename)
        return prompt.format(**kwargs)
    
    def clear_cache(self) -> None:
        """キャッシュをクリア"""
        self._cache.clear()


# シングルトンインスタンス
_prompt_loader: Optional[PromptLoader] = None


def get_prompt_loader() -> PromptLoader:
    """PromptLoaderのシングルトンインスタンスを取得"""
    global _prompt_loader
    if _prompt_loader is None:
        _prompt_loader = PromptLoader()
    return _prompt_loader

