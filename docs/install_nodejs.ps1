# Node.js インストールスクリプト
# このスクリプトはwingetを使用してNode.js LTSをインストールします

Write-Host "=========================================="
Write-Host "Node.js インストールスクリプト"
Write-Host "=========================================="
Write-Host ""

# wingetの確認
$wingetPath = Get-Command winget -ErrorAction SilentlyContinue
if (-not $wingetPath) {
    Write-Host "エラー: wingetが見つかりません" -ForegroundColor Red
    Write-Host "Windows 10/11の最新バージョンが必要です" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "代替方法: https://nodejs.org/ から手動でインストールしてください" -ForegroundColor Yellow
    exit 1
}

Write-Host "wingetが見つかりました: $($wingetPath.Source)" -ForegroundColor Green
Write-Host ""

# 管理者権限の確認
$isAdmin = ([Security.Principal.WindowsPrincipal] [Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)

if (-not $isAdmin) {
    Write-Host "警告: 管理者権限で実行されていません" -ForegroundColor Yellow
    Write-Host "管理者権限で実行することを推奨します" -ForegroundColor Yellow
    Write-Host ""
    $continue = Read-Host "続行しますか？ (Y/N)"
    if ($continue -ne "Y" -and $continue -ne "y") {
        Write-Host "インストールをキャンセルしました" -ForegroundColor Yellow
        exit 0
    }
}

Write-Host "Node.js LTSをインストールします..." -ForegroundColor Cyan
Write-Host ""

try {
    # wingetでNode.js LTSをインストール
    winget install OpenJS.NodeJS.LTS --accept-package-agreements --accept-source-agreements
    
    Write-Host ""
    Write-Host "=========================================="
    Write-Host "インストールが完了しました！" -ForegroundColor Green
    Write-Host "=========================================="
    Write-Host ""
    Write-Host "次のステップ:" -ForegroundColor Cyan
    Write-Host "1. PowerShellを再起動してください" -ForegroundColor Yellow
    Write-Host "2. 以下のコマンドでインストールを確認:" -ForegroundColor Yellow
    Write-Host "   node --version" -ForegroundColor White
    Write-Host "   npm --version" -ForegroundColor White
    Write-Host ""
    Write-Host "3. プロジェクトの依存関係をインストール:" -ForegroundColor Yellow
    Write-Host "   cd web" -ForegroundColor White
    Write-Host "   npm install" -ForegroundColor White
    Write-Host ""
    
} catch {
    Write-Host ""
    Write-Host "エラー: インストールに失敗しました" -ForegroundColor Red
    Write-Host $_.Exception.Message -ForegroundColor Red
    Write-Host ""
    Write-Host "代替方法: https://nodejs.org/ から手動でインストールしてください" -ForegroundColor Yellow
    exit 1
}

