"""Check Phase Tasks - 分析タスク"""
from crewai import Task
from typing import Dict, Any, Optional
from src.crewai.agents.performance_analyst import create_performance_analyst_agent
from src.crewai.agents.external_benchmark_analyst import create_external_benchmark_analyst_agent
from src.crewai.agents.domain_trend_analyst import create_domain_trend_analyst_agent
from src.crewai.state.article_state import ArticleState
import json


def create_performance_analysis_task(article_state: ArticleState, context: Dict[str, Any]) -> Task:
    """
    パフォーマンス分析タスクを作成
    
    Args:
        article_state: Article State
        context: コンテキスト情報
    
    Returns:
        Performance Analysis Task
    """
    analyst = create_performance_analyst_agent()
    
    qiita_url = article_state.qiita_url
    qiita_item_id = article_state.qiita_item_id
    
    return Task(
        description=f"""
        以下の記事のパフォーマンスを分析してください：
        
        記事ID: {article_state.article_id}
        Qiita URL: {qiita_url or "（未投稿）"}
        Qiita Item ID: {qiita_item_id or "（未投稿）"}
        
        分析項目:
        1. いいね数・閲覧数・コメント数の分析
        2. 時系列でのパフォーマンス追跡（投稿後1日、1週間、1ヶ月）
        3. 高パフォーマンス記事との比較
        4. 低パフォーマンス記事の原因分析
        5. 改善のヒント
        
        分析結果はJSON形式で出力してください：
        {{
            "likes_count": いいね数,
            "views_count": 閲覧数,
            "comments_count": コメント数,
            "performance_score": パフォーマンススコア（1-10）,
            "trend": "上昇/下降/横ばい",
            "high_performance_features": ["特徴1", "特徴2", ...],
            "low_performance_reasons": ["原因1", "原因2", ...],
            "improvement_hints": ["改善ヒント1", "改善ヒント2", ...]
        }}
        """,
        agent=analyst,
        expected_output="パフォーマンス分析レポート（JSON形式）。いいね数、閲覧数、コメント数、パフォーマンススコア、トレンド、特徴、原因、改善ヒントを含む。"
    )


def create_external_benchmark_task(article_state: ArticleState, context: Dict[str, Any]) -> Task:
    """
    外部ベンチマーク分析タスクを作成
    
    Args:
        article_state: Article State
        context: コンテキスト情報
    
    Returns:
        External Benchmark Task
    """
    analyst = create_external_benchmark_analyst_agent()
    
    return Task(
        description=f"""
        以下の記事を、同じタグ・トピックの他記事やトレンド記事と比較分析してください：
        
        記事ID: {article_state.article_id}
        タイトル: {article_state.plan.get('title', '') if article_state.plan else ''}
        内容: {article_state.content[:500] if article_state.content else ''}...
        
        比較項目:
        1. 同じタグ・トピックの他記事との比較
        2. トレンド記事との差異分析
        3. 競合記事の特徴抽出
        4. 自記事の強み・弱み
        5. 改善のヒント
        
        分析結果はJSON形式で出力してください：
        {{
            "comparison_articles": ["比較記事1", "比較記事2", ...],
            "differences": {{
                "title": "タイトルの違い",
                "structure": "構成の違い",
                "content_depth": "内容の深さの違い"
            }},
            "competitor_features": ["競合記事の特徴1", "特徴2", ...],
            "strengths": ["自記事の強み1", "強み2", ...],
            "weaknesses": ["自記事の弱み1", "弱み2", ...],
            "improvement_hints": ["改善ヒント1", "改善ヒント2", ...]
        }}
        """,
        agent=analyst,
        expected_output="外部ベンチマーク分析レポート（JSON形式）。比較記事、差異、競合特徴、強み・弱み、改善ヒントを含む。"
    )


def create_domain_trend_task(article_state: ArticleState, context: Dict[str, Any]) -> Task:
    """
    ドメイントレンド分析タスクを作成
    
    Args:
        article_state: Article State
        context: コンテキスト情報
    
    Returns:
        Domain Trend Task
    """
    analyst = create_domain_trend_analyst_agent()
    
    return Task(
        description=f"""
        バイオインフォマティクス領域のトレンドを分析し、次回記事のトピック候補を提案してください：
        
        現在の記事:
        タイトル: {article_state.plan.get('title', '') if article_state.plan else ''}
        トピック: {article_state.topic}
        
        分析項目:
        1. バイオインフォマティクス領域の最新トレンド
        2. 技術トレンドの時系列分析
        3. 読者のニーズの把握
        4. 次回記事のトピック候補提案（3-5個）
        5. ドメイン知識の更新
        
        分析結果はJSON形式で出力してください：
        {{
            "current_trends": ["トレンド1", "トレンド2", ...],
            "trend_analysis": "トレンド分析の詳細",
            "reader_needs": ["読者のニーズ1", "ニーズ2", ...],
            "next_topic_candidates": [
                {{
                    "topic": "トピック1",
                    "reason": "提案理由",
                    "priority": 優先度（1-10）
                }},
                ...
            ],
            "domain_knowledge_updates": ["知識更新1", "更新2", ...]
        }}
        """,
        agent=analyst,
        expected_output="ドメイントレンド分析レポート（JSON形式）。現在のトレンド、トレンド分析、読者のニーズ、次回トピック候補、ドメイン知識更新を含む。"
    )

