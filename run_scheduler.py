"""スケジューラーを常駐実行するスクリプト"""
import sys
import time
import signal
from pathlib import Path
from dotenv import load_dotenv

# .envファイルを読み込む
load_dotenv()

from src.scheduler.article_scheduler import get_scheduler
from src.scheduler.draft_scheduler import DraftScheduler


def signal_handler(sig, frame):
    """シグナルハンドラー（Ctrl+Cで停止）"""
    print("\n[SCHEDULER] スケジューラーを停止しています...")
    scheduler = get_scheduler()
    scheduler.stop()
    sys.exit(0)


def main():
    """メイン関数"""
    print("="*60)
    print("Qiita記事自動投稿スケジューラー")
    print("="*60)
    print()
    
    # ドラフトスケジューラーの状態を表示
    draft_scheduler = DraftScheduler()
    status = draft_scheduler.get_status()
    
    print("[現在の状態]")
    print(f"  総ドラフト数: {status['total_drafts']}")
    print(f"  投稿済み: {status['published_count']}")
    print(f"  現在のインデックス: {status['current_index']}")
    if status['next_draft']:
        print(f"  次回投稿予定: {Path(status['next_draft']).name}")
    if status['last_published_date']:
        print(f"  最終投稿日時: {status['last_published_date']}")
    print()
    
    # スケジューラーを開始
    scheduler = get_scheduler()
    scheduler.start()
    
    print("[SCHEDULER] スケジューラーが開始されました")
    print("[SCHEDULER] 毎週金曜日に自動投稿が実行されます")
    print("[SCHEDULER] Ctrl+C で停止できます")
    print()
    
    # シグナルハンドラーを登録
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    
    # 常駐実行
    try:
        while True:
            time.sleep(60)  # 1分ごとにチェック
    except KeyboardInterrupt:
        signal_handler(None, None)


if __name__ == "__main__":
    main()
