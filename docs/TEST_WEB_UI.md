# Web UI動作確認ガイド

## 前提条件

### 1. 環境変数の設定

プロジェクトルートに `.env` ファイルを作成し、以下を設定してください：

```env
# Google Gemini API キー（必須）
GEMINI_API_KEY=your_gemini_api_key_here

# Qiita API アクセストークン（投稿時のみ必須）
QIITA_ACCESS_TOKEN=your_qiita_access_token_here

# Serper API キー（Web検索機能を使用する場合、オプション）
SERPER_API_KEY=your_serper_api_key_here
```

### 2. 設定ファイルの確認

`config/config.json` が存在することを確認してください。存在しない場合は：

```bash
cp config/config.example.json config/config.json
```

## セットアップ手順

### 1. FastAPIバックエンドのセットアップ

```bash
# 仮想環境の作成（推奨）
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

# 依存関係のインストール
pip install -r requirements.txt
pip install -r api/requirements.txt
```

### 2. Next.jsフロントエンドのセットアップ

```bash
cd web
npm install
```

## 起動方法

### 1. FastAPIサーバーの起動

```bash
# プロジェクトルートから
uvicorn api.app.main:app --reload --host 0.0.0.0 --port 8000
```

または

```bash
cd api
uvicorn app.main:app --reload
```

サーバーは `http://localhost:8000` で起動します。

### 2. Next.js開発サーバーの起動

```bash
cd web
npm run dev
```

サーバーは `http://localhost:3000` で起動します。

## 動作確認

### 1. APIヘルスチェック

ブラウザまたはcurlで以下にアクセス：

```bash
curl http://localhost:8000/health
```

期待されるレスポンス：
```json
{"status": "healthy"}
```

### 2. APIドキュメント確認

ブラウザで以下にアクセス：
- Swagger UI: http://localhost:8000/docs
- ReDoc: http://localhost:8000/redoc

### 3. Web UI確認

ブラウザで以下にアクセス：
- ダッシュボード: http://localhost:3000
- 新規記事作成: http://localhost:3000/new

## テスト手順

### テスト1: 記事一覧取得

1. ブラウザで `http://localhost:3000` にアクセス
2. ダッシュボードに既存の記事が表示されることを確認
3. 記事がない場合は「記事がありません」と表示されることを確認

### テスト2: 新規記事作成

1. 「新規記事作成」ボタンをクリック
2. 記事のテーマを入力（例: "PythonでQiita記事投稿を自動化するパイプライン構築"）
3. 対象読者を選択
4. オプション設定を展開して設定（任意）
5. 「記事を作成」ボタンをクリック
6. Plan Phaseが実行され、記事詳細ページに遷移することを確認

### テスト3: 記事詳細表示

1. ダッシュボードから記事カードをクリック
2. 記事詳細ページが表示されることを確認
3. 以下のタブが表示されることを確認：
   - Article: Markdownエディタとプレビュー
   - AI Process: Phase実行ボタン
   - Feedback: Human Feedback入力パネル
   - Analysis: 分析結果表示

### テスト4: Markdown編集

1. 記事詳細ページの「Article」タブを開く
2. 「編集」ボタンをクリック
3. Markdownを編集
4. 「保存」ボタンをクリック
5. 変更が保存されることを確認

### テスト5: Phase実行

1. 記事詳細ページの「AI Process」タブを開く
2. 「Plan Phase実行」ボタンをクリック
3. 実行中は「実行中...」と表示されることを確認
4. 実行完了後、記事状態が更新されることを確認

### テスト6: Human Feedback追加

1. 記事詳細ページの「Feedback」タブを開く
2. フィードバック内容を入力
3. 対象セクション、意図、優先度を設定（任意）
4. 「フィードバックを送信」ボタンをクリック
5. フィードバック履歴に追加されることを確認

### テスト7: Qiita投稿

1. 記事詳細ページで「Qiitaに投稿」ボタンをクリック
2. 確認ダイアログで「OK」をクリック
3. 投稿完了後、「Qiitaで見る」ボタンが表示されることを確認

## トラブルシューティング

### エラー: CORSエラー

FastAPIのCORS設定を確認してください。`api/app/main.py` で `allow_origins` に `http://localhost:3000` が含まれていることを確認。

### エラー: API接続エラー

1. FastAPIサーバーが起動していることを確認
2. `web/lib/api.ts` の `API_BASE_URL` が正しいことを確認
3. ブラウザの開発者ツールでネットワークエラーを確認

### エラー: 記事が見つかりません

`data/state/articles/` ディレクトリに記事ファイルが存在することを確認。

### エラー: Module not found

依存関係が正しくインストールされているか確認：

```bash
# Python
pip install -r requirements.txt
pip install -r api/requirements.txt

# Node.js
cd web
npm install
```

### エラー: 環境変数が見つかりません

`.env` ファイルがプロジェクトルートに存在し、必要な環境変数が設定されていることを確認。

## パフォーマンステスト

### APIレスポンス時間

```bash
# 記事一覧取得
time curl http://localhost:8000/api/articles

# 記事詳細取得
time curl http://localhost:8000/api/articles/{article_id}
```

### フロントエンドパフォーマンス

ブラウザの開発者ツールで以下を確認：
- Network タブ: API呼び出しのレスポンス時間
- Performance タブ: ページ読み込み時間
- Console タブ: エラーや警告

## 次のステップ

実装が完了したら、以下を検討してください：

1. **本番環境へのデプロイ**
   - FastAPI: Gunicorn + Nginx
   - Next.js: Vercel または Docker

2. **テストの自動化**
   - APIテスト: pytest
   - E2Eテスト: Playwright または Cypress

3. **監視とログ**
   - ログ集約: ELK Stack または CloudWatch
   - エラー追跡: Sentry

4. **セキュリティ**
   - API認証: JWT または OAuth2
   - レート制限: FastAPI-limiter

