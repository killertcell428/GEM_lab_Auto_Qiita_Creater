"""FastAPI動作確認用の簡単なテストスクリプト"""
import requests
import json
import sys
from typing import Dict, Any

API_BASE_URL = "http://localhost:8000"


def test_health():
    """ヘルスチェックテスト"""
    print("=" * 50)
    print("テスト1: ヘルスチェック")
    print("=" * 50)
    try:
        response = requests.get(f"{API_BASE_URL}/health")
        response.raise_for_status()
        print(f"✓ ステータスコード: {response.status_code}")
        print(f"✓ レスポンス: {response.json()}")
        return True
    except Exception as e:
        print(f"✗ エラー: {e}")
        return False


def test_list_articles():
    """記事一覧取得テスト"""
    print("\n" + "=" * 50)
    print("テスト2: 記事一覧取得")
    print("=" * 50)
    try:
        response = requests.get(f"{API_BASE_URL}/api/articles")
        response.raise_for_status()
        articles = response.json()
        print(f"✓ ステータスコード: {response.status_code}")
        print(f"✓ 記事数: {len(articles)}")
        if articles:
            print(f"✓ 最初の記事: {articles[0].get('title', 'N/A')}")
        return True
    except Exception as e:
        print(f"✗ エラー: {e}")
        return False


def test_create_article():
    """記事作成テスト"""
    print("\n" + "=" * 50)
    print("テスト3: 記事作成")
    print("=" * 50)
    try:
        payload = {
            "topic": "テスト記事: FastAPIとNext.jsの統合テスト",
            "target_audience": "中級者",
            "article_tone": "technical",
            "code_ratio": "medium",
            "theory_depth": "medium"
        }
        response = requests.post(
            f"{API_BASE_URL}/api/articles",
            json=payload
        )
        response.raise_for_status()
        article = response.json()
        print(f"✓ ステータスコード: {response.status_code}")
        print(f"✓ 記事ID: {article.get('id')}")
        print(f"✓ タイトル: {article.get('title')}")
        print(f"✓ Phase: {article.get('phase')}")
        return article.get('id')
    except Exception as e:
        print(f"✗ エラー: {e}")
        if hasattr(e, 'response') and e.response is not None:
            try:
                error_detail = e.response.json()
                print(f"  エラー詳細: {error_detail}")
            except:
                print(f"  レスポンス: {e.response.text}")
        return None


def test_get_article(article_id: str):
    """記事詳細取得テスト"""
    print("\n" + "=" * 50)
    print("テスト4: 記事詳細取得")
    print("=" * 50)
    try:
        response = requests.get(f"{API_BASE_URL}/api/articles/{article_id}")
        response.raise_for_status()
        article = response.json()
        print(f"✓ ステータスコード: {response.status_code}")
        print(f"✓ 記事ID: {article.get('id')}")
        print(f"✓ タイトル: {article.get('title')}")
        print(f"✓ Phase: {article.get('phase')}")
        print(f"✓ ステータス: {article.get('uiStatusText')}")
        return True
    except Exception as e:
        print(f"✗ エラー: {e}")
        return False


def test_update_article(article_id: str):
    """記事更新テスト"""
    print("\n" + "=" * 50)
    print("テスト5: 記事更新")
    print("=" * 50)
    try:
        payload = {
            "content": "# テスト記事\n\nこれはテスト用の記事です。\n\n## セクション1\n\nテストコンテンツ。",
            "title": "更新されたテスト記事"
        }
        response = requests.put(
            f"{API_BASE_URL}/api/articles/{article_id}",
            json=payload
        )
        response.raise_for_status()
        article = response.json()
        print(f"✓ ステータスコード: {response.status_code}")
        print(f"✓ タイトル: {article.get('title')}")
        print(f"✓ コンテンツ長: {len(article.get('markdown', ''))} 文字")
        return True
    except Exception as e:
        print(f"✗ エラー: {e}")
        return False


def test_create_feedback(article_id: str):
    """Human Feedback作成テスト"""
    print("\n" + "=" * 50)
    print("テスト6: Human Feedback作成")
    print("=" * 50)
    try:
        payload = {
            "content": "この記事をもっと詳しく説明してください",
            "target_section": "全体",
            "intent": "もっと詳しく",
            "priority": 7
        }
        response = requests.post(
            f"{API_BASE_URL}/api/articles/{article_id}/feedback",
            json=payload
        )
        response.raise_for_status()
        feedback = response.json()
        print(f"✓ ステータスコード: {response.status_code}")
        print(f"✓ フィードバックID: {feedback.get('feedback_id')}")
        print(f"✓ コンテンツ: {feedback.get('content')}")
        return True
    except Exception as e:
        print(f"✗ エラー: {e}")
        return False


def test_get_feedback(article_id: str):
    """Feedback履歴取得テスト"""
    print("\n" + "=" * 50)
    print("テスト7: Feedback履歴取得")
    print("=" * 50)
    try:
        response = requests.get(f"{API_BASE_URL}/api/articles/{article_id}/feedback")
        response.raise_for_status()
        feedbacks = response.json()
        print(f"✓ ステータスコード: {response.status_code}")
        print(f"✓ フィードバック数: {len(feedbacks)}")
        return True
    except Exception as e:
        print(f"✗ エラー: {e}")
        return False


def test_settings():
    """設定取得テスト"""
    print("\n" + "=" * 50)
    print("テスト8: 設定取得")
    print("=" * 50)
    try:
        response = requests.get(f"{API_BASE_URL}/api/settings")
        response.raise_for_status()
        settings = response.json()
        print(f"✓ ステータスコード: {response.status_code}")
        print(f"✓ 設定キー数: {len(settings)}")
        return True
    except Exception as e:
        print(f"✗ エラー: {e}")
        return False


def main():
    """メインテスト実行"""
    print("\n" + "=" * 50)
    print("FastAPI動作確認テスト")
    print("=" * 50)
    print(f"\nAPIベースURL: {API_BASE_URL}")
    print("注意: FastAPIサーバーが起動していることを確認してください\n")

    results = []

    # テスト1: ヘルスチェック
    results.append(("ヘルスチェック", test_health()))

    # テスト2: 記事一覧取得
    results.append(("記事一覧取得", test_list_articles()))

    # テスト3: 記事作成
    article_id = test_create_article()
    results.append(("記事作成", article_id is not None))

    if article_id:
        # テスト4: 記事詳細取得
        results.append(("記事詳細取得", test_get_article(article_id)))

        # テスト5: 記事更新
        results.append(("記事更新", test_update_article(article_id)))

        # テスト6: Feedback作成
        results.append(("Feedback作成", test_create_feedback(article_id)))

        # テスト7: Feedback履歴取得
        results.append(("Feedback履歴取得", test_get_feedback(article_id)))

    # テスト8: 設定取得
    results.append(("設定取得", test_settings()))

    # 結果サマリ
    print("\n" + "=" * 50)
    print("テスト結果サマリ")
    print("=" * 50)
    passed = sum(1 for _, result in results if result)
    total = len(results)
    for test_name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {test_name}")
    print(f"\n合計: {passed}/{total} テストが成功しました")

    if passed == total:
        print("\n🎉 すべてのテストが成功しました！")
        return 0
    else:
        print(f"\n⚠️  {total - passed} 個のテストが失敗しました")
        return 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n\nテストが中断されました")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n予期しないエラーが発生しました: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

