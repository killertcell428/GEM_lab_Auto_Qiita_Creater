"""CrewAI版エントリーポイント - CLIでワークフローを実行"""
import argparse
import sys
from pathlib import Path
from dotenv import load_dotenv

# .envファイルを読み込む
load_dotenv()

from src.crewai.orchestrator import PDCAOrchestrator
from src.crewai.state.article_state import ArticleState
from src.crewai.human_loop import HumanLoop


def main():
    """メイン関数"""
    parser = argparse.ArgumentParser(description="Qiita記事投稿自動化パイプライン（CrewAI版）")
    
    subparsers = parser.add_subparsers(dest="command", help="実行するコマンド")
    
    # フルPDCAサイクル実行コマンド
    pdca_parser = subparsers.add_parser("pdca", help="フルPDCAサイクルを実行")
    pdca_parser.add_argument(
        "--topic",
        type=str,
        required=True,
        help="記事のトピック"
    )
    pdca_parser.add_argument(
        "--draft-file",
        type=str,
        help="ドラフトファイルのパス（オプション）"
    )
    pdca_parser.add_argument(
        "--auto-publish",
        action="store_true",
        help="自動的にQiitaに投稿（承認をスキップ）"
    )
    
    # 個別Phase実行コマンド
    plan_parser = subparsers.add_parser("plan", help="Plan Phaseのみ実行")
    plan_parser.add_argument(
        "--article-id",
        type=str,
        help="既存の記事ID（指定しない場合は新規作成）"
    )
    plan_parser.add_argument(
        "--topic",
        type=str,
        required=True,
        help="記事のトピック"
    )
    
    do_parser = subparsers.add_parser("do", help="Do Phaseのみ実行")
    do_parser.add_argument(
        "--article-id",
        type=str,
        required=True,
        help="記事ID"
    )
    
    check_parser = subparsers.add_parser("check", help="Check Phaseのみ実行")
    check_parser.add_argument(
        "--article-id",
        type=str,
        required=True,
        help="記事ID"
    )
    
    act_parser = subparsers.add_parser("act", help="Act Phaseのみ実行")
    act_parser.add_argument(
        "--article-id",
        type=str,
        required=True,
        help="記事ID"
    )
    
    # Human Feedbackコマンド
    feedback_parser = subparsers.add_parser("feedback", help="Human Feedbackを追加")
    feedback_parser.add_argument(
        "--article-id",
        type=str,
        required=True,
        help="記事ID"
    )
    feedback_parser.add_argument(
        "--phase",
        type=str,
        choices=["plan", "do", "check", "act"],
        required=True,
        help="現在のPhase"
    )
    feedback_parser.add_argument(
        "--content",
        type=str,
        required=True,
        help="フィードバック内容（自然言語）"
    )
    feedback_parser.add_argument(
        "--priority",
        type=int,
        default=5,
        help="優先度（1-10、デフォルト: 5）"
    )
    
    # 状態確認コマンド
    state_parser = subparsers.add_parser("state", help="Article Stateを表示")
    state_parser.add_argument(
        "--article-id",
        type=str,
        required=True,
        help="記事ID"
    )
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    try:
        if args.command == "pdca":
            # フルPDCAサイクル実行
            orchestrator = PDCAOrchestrator()
            result = orchestrator.execute_full_pdca_cycle(
                topic=args.topic,
                draft_file=args.draft_file,
                auto_publish=args.auto_publish
            )
            
            if result.get("success"):
                print("\n[OK] PDCAサイクルが正常に完了しました")
                print(f"記事ID: {result.get('article_id')}")
                if result.get("phases", {}).get("publish", {}).get("success"):
                    print(f"Qiita URL: {result.get('phases', {}).get('publish', {}).get('url')}")
            else:
                print(f"\n[ERROR] エラー: {result.get('error')}")
                sys.exit(1)
        
        elif args.command == "plan":
            # Plan Phaseのみ実行
            from src.crewai.crews.plan_crew import execute_plan_phase
            
            if args.article_id:
                article_state = ArticleState.load(args.article_id)
            else:
                from datetime import datetime
                article_id = f"article_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                article_state = ArticleState(article_id=article_id)
            
            article_state = execute_plan_phase(article_state, args.topic)
            print(f"\n[OK] Plan Phase完了: {article_state.article_id}")
        
        elif args.command == "do":
            # Do Phaseのみ実行
            from src.crewai.crews.do_crew import execute_do_phase
            
            article_state = ArticleState.load(args.article_id)
            article_state = execute_do_phase(article_state)
            print(f"\n[OK] Do Phase完了: {article_state.article_id}")
        
        elif args.command == "check":
            # Check Phaseのみ実行
            from src.crewai.crews.check_crew import execute_check_phase
            
            article_state = ArticleState.load(args.article_id)
            article_state = execute_check_phase(article_state)
            print(f"\n[OK] Check Phase完了: {article_state.article_id}")
        
        elif args.command == "act":
            # Act Phaseのみ実行
            from src.crewai.crews.act_crew import execute_act_phase
            from src.crewai.state.pdca_state import PDCAState
            from datetime import datetime
            
            article_state = ArticleState.load(args.article_id)
            cycle_id = f"cycle_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            pdca_state = PDCAState(cycle_id=cycle_id, article_id=args.article_id)
            
            article_state, pdca_state = execute_act_phase(article_state, pdca_state)
            print(f"\n[OK] Act Phase完了: {article_state.article_id}")
        
        elif args.command == "feedback":
            # Human Feedback追加
            human_loop = HumanLoop()
            feedback = human_loop.add_feedback(
                article_id=args.article_id,
                phase=args.phase,
                content=args.content,
                priority=args.priority
            )
            print(f"\n[OK] フィードバックを追加しました: {feedback.feedback_id}")
            print(f"次回のPhase実行時に自動的に処理されます")
        
        elif args.command == "state":
            # Article State表示
            article_state = ArticleState.load(args.article_id)
            import json
            print("\n" + "="*50)
            print(f"Article State: {args.article_id}")
            print("="*50)
            print(json.dumps(article_state.to_dict(), ensure_ascii=False, indent=2))
        
    except KeyboardInterrupt:
        print("\n\n中断されました")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] エラーが発生しました: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

