"""Qiita投稿モジュール - Qiita API v2を使用して記事を投稿"""
import requests
from typing import Dict, Any, Optional
from src.config_loader import get_config


class QiitaPublisher:
    """Qiitaへの記事投稿を行うクラス"""
    
    def __init__(self):
        """
        初期化
        
        Qiita API v2を使用して記事を投稿するために、アクセストークンが必要です。
        アクセストークンは .env ファイルの QIITA_ACCESS_TOKEN に設定してください。
        参考: https://qiita.com/api/v2/docs
        """
        config = get_config()
        qiita_config = config.get("api", {}).get("qiita", {})
        self.base_url = qiita_config.get("base_url", "https://qiita.com/api/v2")
        self.access_token = qiita_config.get("access_token", "")
        self.timeout = qiita_config.get("timeout", 30)
        
        if not self.access_token:
            raise ValueError(
                "QIITA_ACCESS_TOKENが設定されていません。.envファイルを確認してください。\n"
                "Qiita API アクセストークンは https://qiita.com/settings/applications で取得できます。"
            )
    
    def publish_article(self, article_data: Dict[str, Any], private: bool = False, tweet: bool = False) -> Dict[str, Any]:
        """
        記事をQiitaに投稿
        
        Args:
            article_data: 記事データ（title, content, tagsを含む）
            private: 限定公開にするかどうか
            tweet: ツイートするかどうか
            
        Returns:
            Dict[str, Any]: 投稿結果
                - success: bool - 成功したかどうか
                - item_id: str - QiitaのアイテムID
                - url: str - 記事のURL
                - error: str - エラーメッセージ（失敗時）
        """
        title = article_data.get("title", "")
        content = article_data.get("content", "")
        tags = article_data.get("tags", [])
        
        if not title or not content:
            return {
                "success": False,
                "error": "タイトルまたはコンテンツが空です"
            }
        
        # タグをQiitaの形式に変換（最大5個まで）
        # Qiita APIの仕様: タグは最大5個まで
        qiita_tags = [{"name": tag} for tag in tags[:5]]
        
        # リクエストボディ
        payload = {
            "title": title,
            "body": content,
            "tags": qiita_tags,
            "private": private,
            "tweet": tweet
        }
        
        # ヘッダー（Qiita API v2ではBearerトークン認証を使用）
        headers = {
            "Authorization": f"Bearer {self.access_token}",
            "Content-Type": "application/json"
        }
        
        # デバッグ情報（トークンの最初と最後の数文字のみ表示）
        token_preview = f"{self.access_token[:10]}...{self.access_token[-10:]}" if len(self.access_token) > 20 else "***"
        print(f"[DEBUG] アクセストークン: {token_preview}")
        print(f"[DEBUG] リクエストURL: {self.base_url}/items")
        print(f"[DEBUG] ペイロード: title={title[:50]}, tags={len(tags)}個")
        
        try:
            response = requests.post(
                f"{self.base_url}/items",
                json=payload,
                headers=headers,
                timeout=self.timeout
            )
            
            # デバッグ情報
            print(f"[DEBUG] レスポンスステータス: {response.status_code}")
            if response.status_code != 201:
                print(f"[DEBUG] レスポンス本文: {response.text[:500]}")
            
            response.raise_for_status()
            result = response.json()
            
            # メール通知を送信
            try:
                from src.notifications.email_sender import EmailSender
                email_sender = EmailSender()
                email_sender.send_article_notification(
                    article_title=title,
                    article_url=result.get("url", ""),
                    article_content=content[:1000] if len(content) > 1000 else content  # 最初の1000文字を送信
                )
            except Exception as e:
                print(f"[WARN] メール通知の送信に失敗しました（投稿は成功）: {str(e)}")
            
            return {
                "success": True,
                "item_id": result.get("id", ""),
                "url": result.get("url", ""),
                "raw_response": result
            }
        except requests.exceptions.HTTPError as e:
            error_msg = f"HTTPエラー: {e.response.status_code}"
            try:
                error_detail = e.response.json()
                error_msg += f" - {error_detail}"
                
                # 403エラーの場合、より詳細な情報を提供
                if e.response.status_code == 403:
                    error_msg += "\n\n考えられる原因:"
                    error_msg += "\n- アクセストークンのスコープが不足しています（write_qiita スコープが必要）"
                    error_msg += "\n- アクセストークンが無効または期限切れです"
                    error_msg += "\n- アクセストークンが正しく設定されていません"
            except:
                error_msg += f" - {e.response.text}"
            
            return {
                "success": False,
                "error": error_msg
            }
        except requests.exceptions.RequestException as e:
            return {
                "success": False,
                "error": f"リクエストエラー: {str(e)}"
            }
    
    def update_article(self, item_id: str, article_data: Dict[str, Any], private: bool = False) -> Dict[str, Any]:
        """
        既存の記事を更新
        
        Args:
            item_id: QiitaのアイテムID
            article_data: 記事データ
            private: 限定公開にするかどうか
            
        Returns:
            Dict[str, Any]: 更新結果
        """
        title = article_data.get("title", "")
        content = article_data.get("content", "")
        tags = article_data.get("tags", [])
        
        if not title or not content:
            return {
                "success": False,
                "error": "タイトルまたはコンテンツが空です"
            }
        
        # タグをQiitaの形式に変換（最大5個まで）
        if len(tags) > 5:
            tags = tags[:5]
        qiita_tags = [{"name": tag} for tag in tags]
        
        payload = {
            "title": title,
            "body": content,
            "tags": qiita_tags,
            "private": private
        }
        
        headers = {
            "Authorization": f"Bearer {self.access_token}",
            "Content-Type": "application/json"
        }
        
        try:
            response = requests.patch(
                f"{self.base_url}/items/{item_id}",
                json=payload,
                headers=headers,
                timeout=self.timeout
            )
            
            response.raise_for_status()
            result = response.json()
            
            return {
                "success": True,
                "item_id": result.get("id", ""),
                "url": result.get("url", ""),
                "raw_response": result
            }
        except requests.exceptions.RequestException as e:
            return {
                "success": False,
                "error": f"更新エラー: {str(e)}"
            }

