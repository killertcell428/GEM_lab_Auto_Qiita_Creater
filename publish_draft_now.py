"""ドラフト記事を今すぐ投稿するスクリプト（初回実行用）"""
import sys
import argparse
from pathlib import Path
from dotenv import load_dotenv

# .envファイルを読み込む
load_dotenv()

# プロジェクトルートをパスに追加
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from src.scheduler.draft_scheduler import DraftScheduler


def main():
    """メイン関数"""
    parser = argparse.ArgumentParser(description="ドラフト記事を今すぐ投稿")
    parser.add_argument(
        "--yes", "-y",
        action="store_true",
        help="確認なしで投稿"
    )
    args = parser.parse_args()
    """メイン関数"""
    print("="*60)
    print("ドラフト記事の即座投稿")
    print("="*60)
    print()
    
    # ドラフトスケジューラーを初期化
    scheduler = DraftScheduler()
    
    # 現在の状態を表示
    status = scheduler.get_status()
    print("[現在の状態]")
    print(f"  総ドラフト数: {status['total_drafts']}")
    print(f"  投稿済み: {status['published_count']}")
    print(f"  現在のインデックス: {status['current_index']}")
    if status['next_draft']:
        print(f"  次回投稿予定: {Path(status['next_draft']).name}")
    print()
    
    # 確認
    if status['next_draft']:
        draft_name = Path(status['next_draft']).name
        print(f"投稿する記事: {draft_name}")
        print()
        
        if not args.yes:
            try:
                confirm = input("この記事をQiitaに投稿しますか？ (y/N): ")
                if confirm.lower() != 'y':
                    print("投稿をキャンセルしました。")
                    return
            except EOFError:
                # 対話的でない環境（例: パイプ、リダイレクト）の場合は自動実行
                print("[INFO] 対話的でない環境のため、自動的に投稿を実行します。")
                print("[INFO] 確認なしで実行する場合は --yes オプションを使用してください。")
        
        print()
        print("[INFO] 記事を投稿しています...")
        print()
        
        # 記事を投稿
        result = scheduler.publish_next_draft(auto_publish=True)
        
        if result.get("success"):
            print()
            print("="*60)
            print("[SUCCESS] 投稿成功！")
            print("="*60)
            print(f"タイトル: {result.get('title', 'N/A')}")
            print(f"URL: {result.get('url', 'N/A')}")
            print()
            print("[INFO] 注意: 画像が含まれている場合は、Qiitaで編集して画像URLを設定してください。")
        else:
            print()
            print("="*60)
            print("[ERROR] 投稿失敗")
            print("="*60)
            print(f"エラー: {result.get('error', 'Unknown error')}")
            sys.exit(1)
    else:
        print("[ERROR] 投稿するドラフトが見つかりません。")
        print("`data/drafts/docs/blog/` にドラフトファイルが存在するか確認してください。")
        sys.exit(1)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n中断されました")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] エラーが発生しました: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
