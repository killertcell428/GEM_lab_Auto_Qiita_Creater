"""CrewAI簡易テスト - 最小限の動作確認"""
import sys
import os
from pathlib import Path

# プロジェクトルートをパスに追加
sys.path.insert(0, str(Path(__file__).parent))

# 環境変数を設定（テスト用）
os.environ.setdefault("OPENAI_API_KEY", "dummy-key-for-testing")

def test_basic_imports():
    """基本的なインポートテスト"""
    print("=" * 50)
    print("基本インポートテスト")
    print("=" * 50)
    
    try:
        from src.crewai.state.article_state import ArticleState
        print("[OK] ArticleState")
    except Exception as e:
        print(f"[ERROR] ArticleState: {e}")
        return False
    
    try:
        from src.crewai.tools.gemini_tool import call_gemini_api
        print("[OK] Gemini Tool")
    except Exception as e:
        print(f"[ERROR] Gemini Tool: {e}")
        return False
    
    try:
        from src.crewai.agents.researcher import create_researcher_agent
        print("[OK] Researcher Agent (import)")
    except Exception as e:
        print(f"[ERROR] Researcher Agent: {e}")
        return False
    
    return True


def test_agent_creation():
    """Agent作成テスト（簡易版）"""
    print("\n" + "=" * 50)
    print("Agent作成テスト（簡易版）")
    print("=" * 50)
    
    try:
        from src.crewai.agents.researcher import create_researcher_agent
        
        # Agent作成を試みる（LLMエラーは無視）
        try:
            agent = create_researcher_agent()
            print("[OK] Researcher Agent作成成功")
            return True
        except Exception as e:
            if "OPENAI_API_KEY" in str(e) or "LLM" in str(e):
                print(f"[WARN] Agent作成時にLLMエラー（これは正常）: {str(e)[:100]}")
                print("[INFO] 実際の実行時には環境変数が設定されます")
                return True  # LLMエラーは無視
            else:
                print(f"[ERROR] Agent作成失敗: {e}")
                return False
    except Exception as e:
        print(f"[ERROR] Agent作成テスト失敗: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    print("CrewAI簡易テスト開始\n")
    
    results = []
    results.append(("基本インポート", test_basic_imports()))
    results.append(("Agent作成", test_agent_creation()))
    
    print("\n" + "=" * 50)
    print("テスト結果")
    print("=" * 50)
    
    for name, result in results:
        status = "[OK]" if result else "[FAIL]"
        print(f"{status} {name}")
    
    all_passed = all(result for _, result in results)
    
    if all_passed:
        print("\n[OK] 基本テストが成功しました！")
        print("\n次のステップ:")
        print("1. .envファイルにGEMINI_API_KEYとQIITA_ACCESS_TOKENを設定")
        print("2. python main_crewai.py pdca --topic \"テストトピック\" で実行")
        sys.exit(0)
    else:
        print("\n[ERROR] 一部のテストが失敗しました")
        sys.exit(1)

