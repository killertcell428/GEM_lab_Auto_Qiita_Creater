"""APIエンドポイントのテストスクリプト"""
import requests
import json
from typing import Dict, Any

BASE_URL = "http://localhost:8000"

def test_health():
    """ヘルスチェックテスト"""
    print("=" * 50)
    print("ヘルスチェックテスト")
    print("=" * 50)
    response = requests.get(f"{BASE_URL}/health")
    print(f"Status: {response.status_code}")
    print(f"Response: {response.json()}")
    assert response.status_code == 200
    assert response.json()["status"] == "healthy"
    print("✓ ヘルスチェック成功\n")

def test_root():
    """ルートエンドポイントテスト"""
    print("=" * 50)
    print("ルートエンドポイントテスト")
    print("=" * 50)
    response = requests.get(f"{BASE_URL}/")
    print(f"Status: {response.status_code}")
    print(f"Response: {response.json()}")
    assert response.status_code == 200
    print("✓ ルートエンドポイント成功\n")

def test_list_articles():
    """記事一覧取得テスト"""
    print("=" * 50)
    print("記事一覧取得テスト")
    print("=" * 50)
    response = requests.get(f"{BASE_URL}/api/articles")
    print(f"Status: {response.status_code}")
    data = response.json()
    print(f"記事数: {len(data)}")
    if data:
        print(f"最初の記事: {json.dumps(data[0], ensure_ascii=False, indent=2)}")
    assert response.status_code == 200
    assert isinstance(data, list)
    print("✓ 記事一覧取得成功\n")

def test_get_nonexistent_article():
    """存在しない記事取得テスト（エラーハンドリング）"""
    print("=" * 50)
    print("存在しない記事取得テスト（エラーハンドリング）")
    print("=" * 50)
    response = requests.get(f"{BASE_URL}/api/articles/nonexistent_id_12345")
    print(f"Status: {response.status_code}")
    print(f"Response: {response.json()}")
    assert response.status_code == 404
    print("✓ 404エラーハンドリング成功\n")

def test_create_article_validation_error():
    """記事作成バリデーションエラーテスト"""
    print("=" * 50)
    print("記事作成バリデーションエラーテスト")
    print("=" * 50)
    # 必須フィールド（topic）を欠いたリクエスト
    response = requests.post(f"{BASE_URL}/api/articles", json={})
    print(f"Status: {response.status_code}")
    print(f"Response: {response.json()}")
    assert response.status_code == 422
    print("✓ バリデーションエラーハンドリング成功\n")

def test_update_nonexistent_article():
    """存在しない記事更新テスト（エラーハンドリング）"""
    print("=" * 50)
    print("存在しない記事更新テスト（エラーハンドリング）")
    print("=" * 50)
    response = requests.put(
        f"{BASE_URL}/api/articles/nonexistent_id_12345",
        json={"title": "テストタイトル", "content": "テストコンテンツ"}
    )
    print(f"Status: {response.status_code}")
    print(f"Response: {response.json()}")
    assert response.status_code == 404
    print("✓ 404エラーハンドリング成功\n")

def test_delete_nonexistent_article():
    """存在しない記事削除テスト（エラーハンドリング）"""
    print("=" * 50)
    print("存在しない記事削除テスト（エラーハンドリング）")
    print("=" * 50)
    response = requests.delete(f"{BASE_URL}/api/articles/nonexistent_id_12345")
    print(f"Status: {response.status_code}")
    # 204 No Contentの場合はレスポンスボディがない
    if response.status_code == 404:
        print(f"Response: {response.json()}")
    assert response.status_code in [204, 404]
    print("✓ エラーハンドリング成功\n")

def main():
    """メインテスト実行"""
    print("\n" + "=" * 50)
    print("APIエンドポイントテスト開始")
    print("=" * 50 + "\n")
    
    try:
        test_health()
        test_root()
        test_list_articles()
        test_get_nonexistent_article()
        test_create_article_validation_error()
        test_update_nonexistent_article()
        test_delete_nonexistent_article()
        
        print("=" * 50)
        print("すべてのテストが成功しました！")
        print("=" * 50)
    except AssertionError as e:
        print(f"\n❌ テスト失敗: {e}")
        return 1
    except Exception as e:
        print(f"\n❌ 予期しないエラー: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())

