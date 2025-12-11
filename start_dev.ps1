# Windows PowerShell用開発環境起動スクリプト

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Qiita記事自動化パイプライン 開発環境起動" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""

# 環境チェック
Write-Host "環境チェック中..." -ForegroundColor Yellow

# Python環境チェック
$pythonCmd = Get-Command python -ErrorAction SilentlyContinue
if ($pythonCmd) {
    $pythonVersion = python --version 2>&1
    Write-Host "✓ Python: $pythonVersion" -ForegroundColor Green
} else {
    Write-Host "エラー: Pythonがインストールされていません" -ForegroundColor Red
    exit 1
}

# Node.js環境チェック
$nodeCmd = Get-Command node -ErrorAction SilentlyContinue
if ($nodeCmd) {
    $nodeVersion = node --version 2>&1
    Write-Host "✓ Node.js: $nodeVersion" -ForegroundColor Green
} else {
    Write-Host "エラー: Node.jsがインストールされていません" -ForegroundColor Red
    exit 1
}

# .envファイルチェック
if (-not (Test-Path ".env")) {
    Write-Host "警告: .envファイルが見つかりません" -ForegroundColor Yellow
    Write-Host "環境変数を設定してください（GEMINI_API_KEY等）" -ForegroundColor Yellow
}

# 依存関係チェック
Write-Host ""
Write-Host "依存関係チェック中..." -ForegroundColor Yellow

# Python 3.12環境チェック
$python312Found = $false
$python312Cmd = Get-Command python3.12 -ErrorAction SilentlyContinue
if ($python312Cmd) {
    $python312Found = $true
} else {
    $py312Cmd = Get-Command py -ErrorAction SilentlyContinue
    if ($py312Cmd) {
        $pyVersion = py -3.12 --version 2>&1
        if ($LASTEXITCODE -eq 0) {
            $python312Found = $true
        }
    }
    if (-not $python312Found) {
        Write-Host "警告: Python 3.12が見つかりません" -ForegroundColor Yellow
        Write-Host "Python 3.12をインストールしてください" -ForegroundColor Yellow
    }
}

# uv環境チェック
$useUv = $false
$uvCmd = Get-Command uv -ErrorAction SilentlyContinue
if ($uvCmd) {
    $uvVersion = uv --version 2>&1
    Write-Host "✓ uv: $uvVersion" -ForegroundColor Green
    $useUv = $true
} else {
    Write-Host "警告: uvが見つかりません" -ForegroundColor Yellow
    $uvInstallCmd = 'powershell -ExecutionPolicy Bypass -c "irm https://astral.sh/uv/install.ps1 | iex"'
    Write-Host "uvをインストールしてください: $uvInstallCmd" -ForegroundColor Yellow
    Write-Host "通常のpipを使用します..." -ForegroundColor Yellow
}

# Python依存関係
if (-not (Test-Path ".venv")) {
    Write-Host ""
    Write-Host "Python仮想環境を作成中..." -ForegroundColor Cyan
    if ($useUv) {
        if ($python312Found) {
            if (Get-Command python3.12 -ErrorAction SilentlyContinue) {
                uv venv .venv --python 3.12
            } else {
                uv venv .venv --python "py -3.12"
            }
        } else {
            Write-Host "警告: Python 3.12が見つかりません。デフォルトのPythonを使用します" -ForegroundColor Yellow
            uv venv .venv
        }
    } else {
        python -m venv .venv
    }
}

Write-Host ""
Write-Host "Python依存関係をインストール中..." -ForegroundColor Cyan

# 仮想環境内のPythonを使用して依存関係をインストール
$venvPython = ".\.venv\Scripts\python.exe"
if ($useUv) {
    & uv pip install -q -r requirements.txt
    & uv pip install -q -r api\requirements.txt
} else {
    & $venvPython -m pip install -q -r requirements.txt
    & $venvPython -m pip install -q -r api\requirements.txt
}

# Node.js依存関係
if (-not (Test-Path "web\node_modules")) {
    Write-Host ""
    Write-Host "Node.js依存関係をインストール中..." -ForegroundColor Cyan
    Push-Location web
    npm install
    Pop-Location
}

# サーバー起動
Write-Host ""
Write-Host "サーバーを起動します..." -ForegroundColor Green
Write-Host ""

# FastAPIサーバーをバックグラウンドで起動
Write-Host "FastAPIサーバーを起動中 (http://localhost:8000)..." -ForegroundColor Cyan
$fastApiCommand = ".\.venv\Scripts\Activate.ps1; uvicorn api.app.main:app --reload --host 0.0.0.0 --port 8000"
Start-Process powershell -ArgumentList "-NoExit", "-Command", $fastApiCommand -WindowStyle Normal

# 少し待つ
Start-Sleep -Seconds 2

# Next.jsサーバーをバックグラウンドで起動
Write-Host "Next.jsサーバーを起動中 (http://localhost:3000)..." -ForegroundColor Cyan
Push-Location web
Start-Process powershell -ArgumentList "-NoExit", "-Command", "npm run dev" -WindowStyle Normal
Pop-Location

Write-Host ""
Write-Host "==========================================" -ForegroundColor Green
Write-Host "サーバーが起動しました！" -ForegroundColor Green
Write-Host "==========================================" -ForegroundColor Green
Write-Host ""
Write-Host "FastAPI: http://localhost:8000" -ForegroundColor White
Write-Host "  - API Docs: http://localhost:8000/docs" -ForegroundColor Gray
Write-Host "  - Health: http://localhost:8000/health" -ForegroundColor Gray
Write-Host ""
Write-Host "Next.js: http://localhost:3000" -ForegroundColor White
Write-Host "  - ダッシュボード: http://localhost:3000" -ForegroundColor Gray
Write-Host "  - 新規記事作成: http://localhost:3000/new" -ForegroundColor Gray
Write-Host ""
Write-Host "各サーバーのウィンドウを閉じることで停止できます" -ForegroundColor Yellow
Write-Host ""
