# クイックスタートガイド

## 最短で動作確認する手順

### 1. 環境変数の設定

プロジェクトルートに `.env` ファイルを作成：

```env
GEMINI_API_KEY=your_gemini_api_key_here
```

### 2. 依存関係のインストール

```bash
# Python依存関係
pip install -r requirements.txt
pip install -r api/requirements.txt

# Node.js依存関係
cd web
npm install
cd ..
```

### 3. サーバー起動

#### 方法1: 自動起動スクリプト（推奨）

**Windows:**
```bash
start_dev.bat
```

**Linux/Mac:**
```bash
./start_dev.sh
```

#### 方法2: 手動起動

**ターミナル1 - FastAPI:**
```bash
uvicorn api.app.main:app --reload
```

**ターミナル2 - Next.js:**
```bash
cd web
npm run dev
```

### 4. ブラウザで確認

- **ダッシュボード**: http://localhost:3000
- **APIドキュメント**: http://localhost:8000/docs

### 5. 簡単なテスト

```bash
# APIテスト
python test_api.py
```

## トラブルシューティング

### ポートが既に使用されている

```bash
# Windows
netstat -ano | findstr :8000
netstat -ano | findstr :3000

# Linux/Mac
lsof -i :8000
lsof -i :3000
```

### モジュールが見つからない

```bash
# Python
pip install -r requirements.txt
pip install -r api/requirements.txt

# Node.js
cd web
npm install
```

詳細は `TEST_WEB_UI.md` を参照してください。

