"""GitHub APIを使用して画像をアップロードし、raw URLを取得するモジュール"""
import requests
import base64
from pathlib import Path
from typing import Dict, Any, Optional
import os
from datetime import datetime


class GitHubImageUploader:
    """GitHub APIを使用して画像をアップロードするクラス"""
    
    def __init__(self):
        """
        初期化
        
        環境変数から以下を取得:
        - GITHUB_TOKEN: GitHub Personal Access Token
        - GITHUB_REPO: リポジトリ名（例: username/repo）
        - GITHUB_BRANCH: ブランチ名（デフォルト: main）
        """
        self.token = os.getenv("GITHUB_TOKEN", "")
        self.repo = os.getenv("GITHUB_REPO", "")
        self.branch = os.getenv("GITHUB_BRANCH", "main")
        
        if not self.token:
            print("[WARN] GITHUB_TOKENが設定されていません。画像の自動アップロードはスキップされます。")
        if not self.repo:
            print("[WARN] GITHUB_REPOが設定されていません。画像の自動アップロードはスキップされます。")
    
    def is_configured(self) -> bool:
        """GitHub設定が完了しているか確認"""
        return bool(self.token and self.repo)
    
    def upload_image(self, image_path: Path, remote_path: Optional[str] = None) -> Optional[str]:
        """
        画像をGitHubにアップロードしてraw URLを取得
        
        Args:
            image_path: アップロードする画像ファイルのパス
            remote_path: リモートパス（指定しない場合は自動生成）
            
        Returns:
            raw URL（成功時）、None（失敗時）
        """
        if not self.is_configured():
            print(f"[WARN] GitHub設定が不完全です。画像のアップロードをスキップ: {image_path.name}")
            return None
        
        if not image_path.exists():
            print(f"[WARN] 画像ファイルが見つかりません: {image_path}")
            return None
        
        try:
            # リモートパスを生成（指定がない場合）
            if remote_path is None:
                # 日付ベースのディレクトリ構造
                date_str = datetime.now().strftime("%Y/%m")
                remote_path = f"images/{date_str}/{image_path.name}"
            
            # 画像ファイルを読み込む
            with open(image_path, "rb") as f:
                image_data = f.read()
            
            # Base64エンコード
            encoded_content = base64.b64encode(image_data).decode("utf-8")
            
            # GitHub APIでファイルをアップロード
            url = f"https://api.github.com/repos/{self.repo}/contents/{remote_path}"
            
            headers = {
                "Authorization": f"token {self.token}",
                "Accept": "application/vnd.github.v3+json"
            }
            
            # 既存ファイルをチェック
            check_response = requests.get(url, headers=headers)
            sha = None
            if check_response.status_code == 200:
                # 既存ファイルがある場合はSHAを取得
                sha = check_response.json().get("sha")
            
            # ファイルをアップロード/更新
            payload = {
                "message": f"Upload image: {image_path.name}",
                "content": encoded_content,
                "branch": self.branch
            }
            
            if sha:
                payload["sha"] = sha  # 更新の場合
            
            response = requests.put(url, json=payload, headers=headers)
            response.raise_for_status()
            
            # raw URLを生成
            raw_url = f"https://raw.githubusercontent.com/{self.repo}/{self.branch}/{remote_path}"
            
            print(f"[INFO] 画像をアップロードしました: {image_path.name} -> {raw_url}")
            return raw_url
            
        except requests.exceptions.RequestException as e:
            print(f"[ERROR] 画像のアップロードに失敗しました ({image_path.name}): {str(e)}")
            return None
        except Exception as e:
            print(f"[ERROR] 画像のアップロードエラー ({image_path.name}): {str(e)}")
            import traceback
            traceback.print_exc()
            return None
    
    def upload_multiple_images(self, image_paths: list[Path], base_remote_path: str = "images") -> Dict[str, Optional[str]]:
        """
        複数の画像をアップロード
        
        Args:
            image_paths: アップロードする画像ファイルのパスリスト
            base_remote_path: ベースリモートパス
            
        Returns:
            画像ファイル名 -> raw URL の辞書
        """
        results = {}
        
        for image_path in image_paths:
            if not image_path.exists():
                results[image_path.name] = None
                continue
            
            # リモートパスを生成
            date_str = datetime.now().strftime("%Y/%m")
            remote_path = f"{base_remote_path}/{date_str}/{image_path.name}"
            
            raw_url = self.upload_image(image_path, remote_path)
            results[image_path.name] = raw_url
        
        return results
