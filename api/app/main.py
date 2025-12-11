"""FastAPI Application - CrewAI Web API"""
from fastapi import FastAPI, Request, status
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.exceptions import RequestValidationError
from dotenv import load_dotenv
import sys
import traceback
from pathlib import Path
from typing import Dict, Any

# プロジェクトルートをパスに追加
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# .envファイルを読み込む
load_dotenv(project_root / ".env")

from api.app.routers import articles, feedback, phases, streaming, settings, approval, research, metrics
from src.scheduler import get_scheduler

app = FastAPI(
    title="Qiita記事投稿自動化パイプライン API",
    description="CrewAIベースのPDCAサイクルシステムのWeb API",
    version="1.0.0"
)


@app.on_event("startup")
async def startup_event():
    """アプリケーション起動時の処理"""
    try:
        scheduler = get_scheduler()
        scheduler.start()
        print("[STARTUP] スケジューラーを開始しました")
    except Exception as e:
        print(f"[WARN] スケジューラーの開始に失敗しました: {str(e)}")


@app.on_event("shutdown")
async def shutdown_event():
    """アプリケーション終了時の処理"""
    try:
        scheduler = get_scheduler()
        scheduler.stop()
        print("[SHUTDOWN] スケジューラーを停止しました")
    except Exception as e:
        print(f"[WARN] スケジューラーの停止に失敗しました: {str(e)}")

# CORS設定（Next.jsからのアクセスを許可）
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://127.0.0.1:3000"],  # Next.jsのデフォルトポート
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ルーターの登録
app.include_router(articles.router, prefix="/api/articles", tags=["articles"])
app.include_router(feedback.router, prefix="/api/articles", tags=["feedback"])
app.include_router(phases.router, prefix="/api/articles", tags=["phases"])
app.include_router(streaming.router, prefix="/api/articles", tags=["streaming"])
app.include_router(settings.router, prefix="/api", tags=["settings"])
app.include_router(approval.router, prefix="/api/articles", tags=["approval"])
app.include_router(research.router, prefix="/api/articles", tags=["research"])
app.include_router(metrics.router, prefix="/api/articles", tags=["metrics"])


@app.get("/")
async def root():
    """ルートエンドポイント"""
    return {
        "message": "Qiita記事投稿自動化パイプライン API",
        "version": "1.0.0",
        "docs": "/docs"
    }


@app.get("/health")
async def health():
    """ヘルスチェック"""
    return {"status": "healthy"}


# グローバル例外ハンドラー
@app.exception_handler(Exception)
async def global_exception_handler(request: Request, exc: Exception):
    """グローバル例外ハンドラー"""
    error_detail = {
        "error": type(exc).__name__,
        "message": str(exc),
        "path": str(request.url.path)
    }
    
    # デバッグモードではスタックトレースも含める
    import os
    if os.getenv("DEBUG", "false").lower() == "true":
        error_detail["traceback"] = traceback.format_exc()
    
    print(f"[ERROR] グローバル例外ハンドラー: {error_detail}")
    
    return JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content=error_detail
    )


@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    """バリデーションエラーハンドラー"""
    error_detail = {
        "error": "ValidationError",
        "message": "リクエストのバリデーションに失敗しました",
        "details": exc.errors(),
        "path": str(request.url.path)
    }
    
    print(f"[ERROR] バリデーションエラー: {error_detail}")
    
    return JSONResponse(
        status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
        content=error_detail
    )

