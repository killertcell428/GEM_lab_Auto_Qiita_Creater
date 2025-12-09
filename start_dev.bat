@echo off
REM Windows用開発環境起動スクリプト

echo ==========================================
echo Qiita記事自動化パイプライン 開発環境起動
echo ==========================================

REM 環境チェック
echo 環境チェック中...

REM Python環境チェック
python --version >nul 2>&1
if errorlevel 1 (
    echo エラー: Pythonがインストールされていません
    exit /b 1
)

REM Node.js環境チェック
node --version >nul 2>&1
if errorlevel 1 (
    echo エラー: Node.jsがインストールされていません
    exit /b 1
)

REM .envファイルチェック
if not exist ".env" (
    echo 警告: .envファイルが見つかりません
    echo 環境変数を設定してください（GEMINI_API_KEY等）
)

REM 依存関係チェック
echo 依存関係チェック中...

REM Python 3.12環境チェック
python3.12 --version >nul 2>&1
if errorlevel 1 (
    py -3.12 --version >nul 2>&1
    if errorlevel 1 (
        echo 警告: Python 3.12が見つかりません
        echo Python 3.12をインストールしてください
    )
)

REM uv環境チェック
uv --version >nul 2>&1
if errorlevel 1 (
    echo 警告: uvが見つかりません
    echo uvをインストールしてください: powershell -ExecutionPolicy Bypass -c "irm https://astral.sh/uv/install.ps1 | iex"
    echo 通常のpipを使用します...
    set USE_UV=0
) else (
    set USE_UV=1
)

REM Python依存関係
if not exist ".venv" (
    echo Python仮想環境を作成中...
    if "%USE_UV%"=="1" (
        python3.12 --version >nul 2>&1
        if not errorlevel 1 (
            uv venv .venv --python 3.12
        ) else (
            py -3.12 --version >nul 2>&1
            if not errorlevel 1 (
                uv venv .venv --python "py -3.12"
            ) else (
                echo 警告: Python 3.12が見つかりません。デフォルトのPythonを使用します
                uv venv .venv
            )
        )
    ) else (
        python -m venv .venv
    )
)

echo Python依存関係をインストール中...
call .venv\Scripts\activate.bat
if "%USE_UV%"=="1" (
    uv pip install -q -r requirements.txt
    uv pip install -q -r api\requirements.txt
) else (
    pip install -q -r requirements.txt
    pip install -q -r api\requirements.txt
)

REM Node.js依存関係
if not exist "web\node_modules" (
    echo Node.js依存関係をインストール中...
    cd web
    call npm install
    cd ..
)

REM サーバー起動
echo.
echo サーバーを起動します...
echo.

REM FastAPIサーバーをバックグラウンドで起動
echo FastAPIサーバーを起動中 (http://localhost:8000)...
start "FastAPI Server" cmd /k ".venv\Scripts\activate.bat && uvicorn api.app.main:app --reload --host 0.0.0.0 --port 8000"

REM 少し待つ
timeout /t 2 /nobreak >nul

REM Next.jsサーバーをバックグラウンドで起動
echo Next.jsサーバーを起動中 (http://localhost:3000)...
cd web
start "Next.js Server" cmd /k "npm run dev"
cd ..

echo.
echo ==========================================
echo サーバーが起動しました！
echo ==========================================
echo.
echo FastAPI: http://localhost:8000
echo   - API Docs: http://localhost:8000/docs
echo   - Health: http://localhost:8000/health
echo.
echo Next.js: http://localhost:3000
echo   - ダッシュボード: http://localhost:3000
echo   - 新規記事作成: http://localhost:3000/new
echo.
echo 各サーバーのウィンドウを閉じることで停止できます
echo.
pause

