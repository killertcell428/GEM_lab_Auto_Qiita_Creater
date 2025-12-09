# uvセットアップスクリプト
# Python 3.12とuvを使用してプロジェクトをセットアップします

Write-Host "=========================================="
Write-Host "uvセットアップスクリプト"
Write-Host "=========================================="
Write-Host ""

# Python 3.12の確認
Write-Host "Python 3.12を確認中..." -ForegroundColor Cyan
$python312 = Get-Command python3.12 -ErrorAction SilentlyContinue
if (-not $python312) {
    $python312 = Get-Command py -ErrorAction SilentlyContinue
    if ($python312) {
        $version = & py -3.12 --version 2>&1
        if ($LASTEXITCODE -ne 0) {
            Write-Host "エラー: Python 3.12が見つかりません" -ForegroundColor Red
            Write-Host ""
            Write-Host "Python 3.12をインストールしてください:" -ForegroundColor Yellow
            Write-Host "  1. https://www.python.org/downloads/ にアクセス" -ForegroundColor White
            Write-Host "  2. Python 3.12.xをダウンロードしてインストール" -ForegroundColor White
            Write-Host "  3. インストール時に「Add Python to PATH」にチェック" -ForegroundColor White
            Write-Host ""
            Write-Host "または winget を使用:" -ForegroundColor Yellow
            Write-Host "  winget install Python.Python.3.12" -ForegroundColor White
            exit 1
        }
    } else {
        Write-Host "エラー: Python 3.12が見つかりません" -ForegroundColor Red
        exit 1
    }
}

Write-Host "✓ Python 3.12が見つかりました" -ForegroundColor Green
Write-Host ""

# uvの確認
Write-Host "uvを確認中..." -ForegroundColor Cyan
$uv = Get-Command uv -ErrorAction SilentlyContinue
if (-not $uv) {
    Write-Host "uvが見つかりません。インストールします..." -ForegroundColor Yellow
    Write-Host ""
    
    try {
        # uvのインストール
        powershell -ExecutionPolicy Bypass -c "irm https://astral.sh/uv/install.ps1 | iex"
        
        # PATHを更新
        $env:PATH = [System.Environment]::GetEnvironmentVariable("Path","Machine") + ";" + [System.Environment]::GetEnvironmentVariable("Path","User")
        
        # 再確認
        $uv = Get-Command uv -ErrorAction SilentlyContinue
        if (-not $uv) {
            Write-Host "警告: uvのインストール後、PowerShellを再起動してください" -ForegroundColor Yellow
            Write-Host "または、手動で以下を実行:" -ForegroundColor Yellow
            Write-Host "  powershell -ExecutionPolicy Bypass -c `"irm https://astral.sh/uv/install.ps1 | iex`"" -ForegroundColor White
            exit 1
        }
    } catch {
        Write-Host "エラー: uvのインストールに失敗しました" -ForegroundColor Red
        Write-Host $_.Exception.Message -ForegroundColor Red
        exit 1
    }
}

Write-Host "✓ uvが見つかりました: $($uv.Source)" -ForegroundColor Green
Write-Host ""

# 仮想環境の作成
Write-Host "仮想環境を作成中..." -ForegroundColor Cyan
if (Test-Path .venv) {
    Write-Host "警告: .venvディレクトリが既に存在します" -ForegroundColor Yellow
    $overwrite = Read-Host "上書きしますか？ (Y/N)"
    if ($overwrite -eq "Y" -or $overwrite -eq "y") {
        Remove-Item -Recurse -Force .venv
    } else {
        Write-Host "既存の仮想環境を使用します" -ForegroundColor Yellow
    }
}

if (-not (Test-Path .venv)) {
    try {
        # Python 3.12で仮想環境を作成
        if (Get-Command python3.12 -ErrorAction SilentlyContinue) {
            uv venv .venv --python 3.12
        } else {
            uv venv .venv --python "py -3.12"
        }
        Write-Host "✓ 仮想環境を作成しました" -ForegroundColor Green
    } catch {
        Write-Host "エラー: 仮想環境の作成に失敗しました" -ForegroundColor Red
        Write-Host $_.Exception.Message -ForegroundColor Red
        exit 1
    }
}
Write-Host ""

# 依存関係のインストール
Write-Host "依存関係をインストール中..." -ForegroundColor Cyan
try {
    # 仮想環境を有効化して依存関係をインストール
    & .\.venv\Scripts\python.exe -m pip install --upgrade pip
    
    # requirements.txtからインストール
    & .\.venv\Scripts\python.exe -m pip install -r requirements.txt
    & .\.venv\Scripts\python.exe -m pip install -r api\requirements.txt
    
    Write-Host "✓ 依存関係のインストールが完了しました" -ForegroundColor Green
} catch {
    Write-Host "エラー: 依存関係のインストールに失敗しました" -ForegroundColor Red
    Write-Host $_.Exception.Message -ForegroundColor Red
    Write-Host ""
    Write-Host "手動でインストールしてください:" -ForegroundColor Yellow
    Write-Host "  .\.venv\Scripts\Activate.ps1" -ForegroundColor White
    Write-Host "  uv pip install -r requirements.txt" -ForegroundColor White
    exit 1
}

Write-Host ""
Write-Host "=========================================="
Write-Host "セットアップが完了しました！" -ForegroundColor Green
Write-Host "=========================================="
Write-Host ""
Write-Host "次のステップ:" -ForegroundColor Cyan
Write-Host "1. 仮想環境を有効化:" -ForegroundColor Yellow
Write-Host "   .\.venv\Scripts\Activate.ps1" -ForegroundColor White
Write-Host ""
Write-Host "2. サーバーを起動:" -ForegroundColor Yellow
Write-Host "   uvicorn api.app.main:app --reload --host 0.0.0.0 --port 8000" -ForegroundColor White
Write-Host ""
Write-Host "または、start_dev.bat を実行してください" -ForegroundColor Yellow
Write-Host ""

