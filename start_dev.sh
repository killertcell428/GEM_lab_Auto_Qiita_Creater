#!/bin/bash
# 開発環境起動スクリプト

echo "=========================================="
echo "Qiita記事自動化パイプライン 開発環境起動"
echo "=========================================="

# カラー定義
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# 環境チェック
echo -e "${YELLOW}環境チェック中...${NC}"

# Python環境チェック
if ! command -v python &> /dev/null; then
    echo "エラー: Pythonがインストールされていません"
    exit 1
fi

# Node.js環境チェック
if ! command -v node &> /dev/null; then
    echo "エラー: Node.jsがインストールされていません"
    exit 1
fi

# .envファイルチェック
if [ ! -f ".env" ]; then
    echo -e "${YELLOW}警告: .envファイルが見つかりません${NC}"
    echo "環境変数を設定してください（GEMINI_API_KEY等）"
fi

# 依存関係チェック
echo -e "${YELLOW}依存関係チェック中...${NC}"

# Python依存関係
if [ ! -d "venv" ]; then
    echo "Python仮想環境を作成中..."
    python -m venv venv
fi

echo "Python依存関係をインストール中..."
source venv/bin/activate
pip install -q -r requirements.txt
pip install -q -r api/requirements.txt

# Node.js依存関係
if [ ! -d "web/node_modules" ]; then
    echo "Node.js依存関係をインストール中..."
    cd web
    npm install
    cd ..
fi

# サーバー起動
echo -e "${GREEN}サーバーを起動します...${NC}"
echo ""

# FastAPIサーバーをバックグラウンドで起動
echo "FastAPIサーバーを起動中 (http://localhost:8000)..."
source venv/bin/activate
uvicorn api.app.main:app --reload --host 0.0.0.0 --port 8000 > /tmp/fastapi.log 2>&1 &
FASTAPI_PID=$!

# 少し待つ
sleep 2

# Next.jsサーバーをバックグラウンドで起動
echo "Next.jsサーバーを起動中 (http://localhost:3000)..."
cd web
npm run dev > /tmp/nextjs.log 2>&1 &
NEXTJS_PID=$!
cd ..

echo ""
echo -e "${GREEN}=========================================="
echo "サーバーが起動しました！"
echo "==========================================${NC}"
echo ""
echo "FastAPI: http://localhost:8000"
echo "  - API Docs: http://localhost:8000/docs"
echo "  - Health: http://localhost:8000/health"
echo ""
echo "Next.js: http://localhost:3000"
echo "  - ダッシュボード: http://localhost:3000"
echo "  - 新規記事作成: http://localhost:3000/new"
echo ""
echo "ログ:"
echo "  - FastAPI: /tmp/fastapi.log"
echo "  - Next.js: /tmp/nextjs.log"
echo ""
echo "停止するには Ctrl+C を押してください"
echo ""

# クリーンアップ関数
cleanup() {
    echo ""
    echo "サーバーを停止中..."
    kill $FASTAPI_PID 2>/dev/null
    kill $NEXTJS_PID 2>/dev/null
    echo "サーバーを停止しました"
    exit 0
}

# シグナルハンドラ
trap cleanup SIGINT SIGTERM

# 待機
wait

