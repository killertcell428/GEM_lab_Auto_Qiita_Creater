"""メトリクスAPIルーター"""
from fastapi import APIRouter, HTTPException
from typing import Dict, Any, List, Optional
from datetime import datetime, timedelta
from pathlib import Path
from src.crewai.state.article_state import ArticleState

router = APIRouter()


@router.get("/{article_id}/metrics")
async def get_article_metrics(article_id: str) -> Dict[str, Any]:
    """記事のパフォーマンスメトリクス取得（時系列データ含む）"""
    try:
        from src.config_loader import get_config
        from src.publisher.qiita_fetcher import QiitaFetcher
        
        config = get_config()
        crewai_config = config.get("crewai", {})
        state_config = crewai_config.get("state", {})
        articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
        file_path = articles_dir / f"{article_id}.json"
        
        if not file_path.exists():
            raise HTTPException(status_code=404, detail="記事が見つかりません")
        
        article_state = ArticleState.load(article_id)
        
        # 基本メトリクス
        metrics = {
            "article_id": article_id,
            "kpi": article_state.kpi or {},
            "analysis_results": article_state.analysis_results or {},
            "qiita_url": article_state.qiita_url,
            "qiita_item_id": article_state.qiita_item_id
        }
        
        # Qiitaに投稿済みの場合、最新のメトリクスを取得
        if article_state.qiita_item_id:
            try:
                fetcher = QiitaFetcher()
                qiita_item = fetcher.get_item(article_state.qiita_item_id)
                
                if qiita_item:
                    metrics["current_metrics"] = {
                        "likes_count": qiita_item.get("likes_count", 0),
                        "page_views_count": qiita_item.get("page_views_count", 0),
                        "comments_count": qiita_item.get("comments_count", 0),
                        "stocks_count": qiita_item.get("stocks_count", 0),
                        "updated_at": qiita_item.get("updated_at")
                    }
            except Exception as e:
                print(f"[WARN] Qiitaメトリクス取得エラー: {str(e)}")
        
        # 時系列データ（将来の拡張用）
        metrics["time_series"] = []
        
        return metrics
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"メトリクスの取得に失敗しました: {str(e)}")


@router.get("/dashboard/summary")
async def get_dashboard_summary() -> Dict[str, Any]:
    """ダッシュボード用サマリー取得"""
    try:
        from src.config_loader import get_config
        
        config = get_config()
        crewai_config = config.get("crewai", {})
        state_config = crewai_config.get("state", {})
        articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
        
        if not articles_dir.exists():
            return {
                "total_articles": 0,
                "total_likes": 0,
                "total_views": 0,
                "avg_engagement_rate": 0.0,
                "pending_approval_count": 0
            }
        
        total_likes = 0
        total_views = 0
        total_articles = 0
        pending_approval_count = 0
        
        for json_file in articles_dir.glob("*.json"):
            try:
                article_state = ArticleState.load(json_file.stem)
                
                if article_state.pending_approval and article_state.approval_status == "pending":
                    pending_approval_count += 1
                
                if article_state.qiita_item_id:
                    total_articles += 1
                    kpi = article_state.kpi or {}
                    total_likes += kpi.get("likes_count", 0)
                    total_views += kpi.get("page_views_count", 0)
            except Exception as e:
                print(f"Warning: Failed to load article {json_file.stem}: {e}")
                continue
        
        avg_engagement_rate = (total_likes / total_views * 100) if total_views > 0 else 0.0
        
        return {
            "total_articles": total_articles,
            "total_likes": total_likes,
            "total_views": total_views,
            "avg_engagement_rate": round(avg_engagement_rate, 2),
            "pending_approval_count": pending_approval_count
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"ダッシュボードサマリーの取得に失敗しました: {str(e)}")


@router.get("/{article_id}/keywords")
async def get_article_keywords(article_id: str) -> Dict[str, Any]:
    """記事のキーワード抽出と分析"""
    try:
        from src.config_loader import get_config
        import re
        from collections import Counter
        
        config = get_config()
        crewai_config = config.get("crewai", {})
        state_config = crewai_config.get("state", {})
        articles_dir = Path(state_config.get("articles_dir", "data/state/articles"))
        file_path = articles_dir / f"{article_id}.json"
        
        if not file_path.exists():
            raise HTTPException(status_code=404, detail="記事が見つかりません")
        
        article_state = ArticleState.load(article_id)
        
        # キーワード抽出（簡易版）
        keywords = []
        if article_state.content:
            # 日本語の単語を抽出（簡易版）
            text = article_state.content
            # コードブロックを除去
            text = re.sub(r'```[\s\S]*?```', '', text)
            text = re.sub(r'`[^`]+`', '', text)
            # マークダウン記号を除去
            text = re.sub(r'[#*_\[\]()]', '', text)
            # 単語を抽出（2文字以上の日本語）
            words = re.findall(r'[\u3040-\u309F\u30A0-\u30FF\u4E00-\u9FAF]{2,}', text)
            # 頻度をカウント
            word_counts = Counter(words)
            # 上位20個を取得
            keywords = [{"text": word, "weight": count} for word, count in word_counts.most_common(20)]
        
        # 分析結果からもキーワードを抽出
        analysis_keywords = []
        if article_state.analysis_results:
            analysis_str = str(article_state.analysis_results)
            words = re.findall(r'[\u3040-\u309F\u30A0-\u30FF\u4E00-\u9FAF]{2,}', analysis_str)
            word_counts = Counter(words)
            analysis_keywords = [{"text": word, "weight": count} for word, count in word_counts.most_common(10)]
        
        return {
            "article_id": article_id,
            "content_keywords": keywords,
            "analysis_keywords": analysis_keywords,
            "combined_keywords": keywords + analysis_keywords
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"キーワード抽出に失敗しました: {str(e)}")

