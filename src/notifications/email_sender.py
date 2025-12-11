"""メール通知モジュール - 記事投稿時にメールを送信"""
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from typing import List, Optional
from src.config_loader import get_config


class EmailSender:
    """メール送信クラス"""
    
    def __init__(self):
        """初期化"""
        import os
        self.config = get_config()
        email_config = self.config.get("notifications", {}).get("email", {})
        self.enabled = email_config.get("enabled", False)
        self.smtp_server = email_config.get("smtp_server", "smtp.gmail.com")
        self.smtp_port = email_config.get("smtp_port", 587)
        # 環境変数から読み込む（優先）
        self.sender_email = os.getenv("EMAIL_SENDER_ADDRESS", email_config.get("sender_email", ""))
        self.sender_password = os.getenv("EMAIL_SENDER_PASSWORD", email_config.get("sender_password", ""))
        self.recipients = email_config.get("recipients", [])
    
    def send_article_notification(self, article_title: str, article_url: str, article_content: Optional[str] = None) -> bool:
        """
        記事投稿通知メールを送信
        
        Args:
            article_title: 記事のタイトル
            article_url: 記事のURL
            article_content: 記事の内容（オプション、抜粋を送信）
            
        Returns:
            bool: 送信成功した場合True
        """
        if not self.enabled:
            print("[INFO] メール通知が無効です")
            return False
        
        if not self.recipients:
            print("[WARN] メール送信先が設定されていません")
            return False
        
        if not self.sender_email or not self.sender_password:
            print("[WARN] メール送信者の設定が不完全です")
            return False
        
        try:
            # メール本文を作成
            content_preview = ""
            if article_content:
                # 記事の最初の500文字を抜粋
                content_preview = article_content[:500] + "..." if len(article_content) > 500 else article_content
            
            body = f"""GEM Lab（遺伝生態モンスター研究所）から新しい記事が投稿されました。

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

【タイトル】
{article_title}

【URL】
{article_url}

【内容抜粋】
{content_preview if content_preview else "（内容抜粋なし）"}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

記事を確認する: {article_url}
"""
            
            # メールを作成
            msg = MIMEMultipart()
            msg["From"] = self.sender_email
            msg["To"] = ", ".join(self.recipients)
            msg["Subject"] = f"【GEM Lab】新しい記事が投稿されました: {article_title}"
            
            msg.attach(MIMEText(body, "plain", "utf-8"))
            
            # SMTPサーバーに接続して送信
            with smtplib.SMTP(self.smtp_server, self.smtp_port) as server:
                server.starttls()
                server.login(self.sender_email, self.sender_password)
                server.send_message(msg)
            
            print(f"[OK] メール通知を送信しました: {len(self.recipients)}件")
            return True
            
        except Exception as e:
            print(f"[ERROR] メール送信に失敗しました: {str(e)}")
            return False
    
    def send_approval_notification(self, article_title: str, article_url: str, approval_deadline: Optional[str] = None, article_content: Optional[str] = None) -> bool:
        """
        承認待ち通知メールを送信
        
        Args:
            article_title: 記事のタイトル
            article_url: 記事のURL（Web UI）
            approval_deadline: 承認期限（ISO形式）
            article_content: 記事の内容（オプション、抜粋を送信）
            
        Returns:
            bool: 送信成功した場合True
        """
        if not self.enabled:
            print("[INFO] メール通知が無効です")
            return False
        
        if not self.recipients:
            print("[WARN] メール送信先が設定されていません")
            return False
        
        if not self.sender_email or not self.sender_password:
            print("[WARN] メール送信者の設定が不完全です")
            return False
        
        try:
            # メール本文を作成
            content_preview = ""
            if article_content:
                # 記事の最初の500文字を抜粋
                content_preview = article_content[:500] + "..." if len(article_content) > 500 else article_content
            
            # 承認期限の表示
            deadline_text = "24時間以内"
            if approval_deadline:
                try:
                    from datetime import datetime
                    deadline_dt = datetime.fromisoformat(approval_deadline)
                    deadline_text = deadline_dt.strftime("%Y年%m月%d日 %H:%M")
                except:
                    deadline_text = approval_deadline
            
            body = f"""新しい記事が作成されました。承認をお願いします。

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

【タイトル】
{article_title}

【承認期限】
{deadline_text}

【内容抜粋】
{content_preview if content_preview else "（内容準備中）"}

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

記事を確認・承認する: {article_url}

承認期限を過ぎると自動的に投稿されます。
"""
            
            # メールを作成
            msg = MIMEMultipart()
            msg["From"] = self.sender_email
            msg["To"] = ", ".join(self.recipients)
            msg["Subject"] = f"【承認待ち】{article_title}"
            
            msg.attach(MIMEText(body, "plain", "utf-8"))
            
            # SMTPサーバーに接続して送信
            with smtplib.SMTP(self.smtp_server, self.smtp_port) as server:
                server.starttls()
                server.login(self.sender_email, self.sender_password)
                server.send_message(msg)
            
            print(f"[OK] 承認通知メールを送信しました: {len(self.recipients)}件")
            return True
            
        except Exception as e:
            print(f"[ERROR] メール送信に失敗しました: {str(e)}")
            return False

