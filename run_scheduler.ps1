# Qiita記事自動投稿スケジューラー起動スクリプト（PowerShell版）

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "Qiita記事自動投稿スケジューラー" -ForegroundColor Cyan
Write-Host "============================================================" -ForegroundColor Cyan
Write-Host ""

# 仮想環境を有効化
$venvPath = ".\.venv\Scripts\Activate.ps1"
if (Test-Path $venvPath) {
    Write-Host "[INFO] 仮想環境を有効化しています..." -ForegroundColor Yellow
    & $venvPath
} else {
    Write-Host "[WARN] 仮想環境が見つかりません: $venvPath" -ForegroundColor Yellow
    Write-Host "[WARN] 仮想環境を作成してください: uv venv .venv --python 3.12" -ForegroundColor Yellow
}

Write-Host ""
Write-Host "[INFO] スケジューラーを起動しています..." -ForegroundColor Green
Write-Host "[INFO] 毎週金曜日に自動投稿が実行されます" -ForegroundColor Green
Write-Host "[INFO] Ctrl+C で停止できます" -ForegroundColor Green
Write-Host ""

# Pythonスクリプトを実行
python run_scheduler.py
