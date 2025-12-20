# ドラフト記事を今すぐ投稿するスクリプト（初回実行用）

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host "ドラフト記事の即座投稿" -ForegroundColor Cyan
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
Write-Host "[INFO] 記事を投稿しています..." -ForegroundColor Green
Write-Host ""

# Pythonスクリプトを実行
python publish_draft_now.py
