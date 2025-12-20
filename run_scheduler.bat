@echo off
REM Qiita記事自動投稿スケジューラー起動スクリプト（バッチ版）

echo ============================================================
echo Qiita記事自動投稿スケジューラー
echo ============================================================
echo.

REM 仮想環境を有効化
if exist ".venv\Scripts\activate.bat" (
    echo [INFO] 仮想環境を有効化しています...
    call .venv\Scripts\activate.bat
) else (
    echo [WARN] 仮想環境が見つかりません
    echo [WARN] 仮想環境を作成してください: uv venv .venv --python 3.12
)

echo.
echo [INFO] スケジューラーを起動しています...
echo [INFO] 毎週金曜日に自動投稿が実行されます
echo [INFO] Ctrl+C で停止できます
echo.

REM Pythonスクリプトを実行
python run_scheduler.py

pause
