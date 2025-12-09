"""Do Phase Tasks - 執筆・レビュータスク"""
from crewai import Task
from typing import Dict, Any
from src.crewai.agents.writer import create_writer_agent
from src.crewai.agents.reviewer import create_reviewer_agent
from src.crewai.agents.implementation_verifier import create_implementation_verifier_agent
from src.crewai.agents.knowledge_extractor import create_knowledge_extractor_agent
from src.prompt_loader import get_prompt_loader
import json


def create_writing_task(plan: Dict[str, Any], research_report: str, context: Dict[str, Any]) -> Task:
    """
    執筆タスクを作成
    
    Args:
        plan: 記事構成プラン
        research_report: 調査レポート
        context: コンテキスト情報
    
    Returns:
        Writing Task
    """
    writer = create_writer_agent()
    
    # プロンプトを読み込む
    prompt_loader = get_prompt_loader()
    article_prompt_template = prompt_loader.load_prompt_or_default(
        "article_generation_prompt.md",
        """以下の構成プランに基づいて、Qiita記事を執筆してください。

## 構成プラン
{plan}

## 調査結果
{research_report}

## 過去の投稿
{past_articles}

記事をMarkdown形式で出力してください。
"""
    )
    
    past_articles = context.get("past_articles", "なし")
    trending_features = context.get("trending_features", "")
    
    # ドラフト内容を取得（contextから）
    draft_content = context.get("draft_content", "")
    additional_requirements = context.get("additional_requirements", "")
    
    plan_str = json.dumps(plan, ensure_ascii=False, indent=2)
    
    article_prompt = article_prompt_template.format(
        draft_content=draft_content if draft_content else "（ドラフトなし）",
        additional_requirements=additional_requirements if additional_requirements else "",
        past_articles_summary=past_articles,
        trending_article_features=trending_features
    )
    
    return Task(
        description=f"""
        以下の構成プランに基づいて、Qiita記事を執筆してください：
        
        {article_prompt}
        
        構成プラン（JSON）:
        {plan_str}
        
        調査結果:
        {research_report}
        
        重要: 
        - 記事の最初の行に、タイトルを`# `で始まるMarkdownの見出し形式で記述してください
        - すべてのセクション（導入、概要、技術詳細、実践編、応用、まとめ、FAQ、参考文献）を必ず含めてください
        - 最低3000文字以上、理想的には5000文字以上の充実した内容にしてください
        - コード例は完全に記述し、コピペで動くレベルにしてください
        - 技術解説が主軸であり、世界観の説明にならないよう注意してください
        """,
        agent=writer,
        expected_output="記事本文（Markdown形式）。タイトル、すべてのセクション、コード例を含む完全な記事。"
    )


def create_review_task(article_content: str, plan: Dict[str, Any], context: Dict[str, Any]) -> Task:
    """
    レビュータスクを作成
    
    Args:
        article_content: 記事本文（空文字列の場合は、前のタスクの出力を参照）
        plan: 記事構成プラン
        context: コンテキスト情報
    
    Returns:
        Review Task
    """
    reviewer = create_reviewer_agent()
    
    plan_str = json.dumps(plan, ensure_ascii=False, indent=2)
    
    # article_contentが空の場合は、前のタスク（writing_task）の出力を参照することを明示
    if not article_content:
        article_content_instruction = """
        前のタスク（記事執筆タスク）で生成された記事本文をレビューしてください。
        前のタスクの出力を参照して、その記事の品質をチェックしてください。
        """
    else:
        article_content_instruction = f"""
        以下の記事をレビューしてください：
        
        記事本文:
        {article_content[:2000]}...（全文は前のタスクの出力を参照）
        """
    
    return Task(
        description=f"""
        {article_content_instruction}
        
        構成プラン:
        {plan_str}
        
        チェック項目:
        1. 技術的な内容の正確性
        2. 構成の論理性（導入→説明→実装→まとめの流れ）
        3. コード例の完全性（コピペで動くレベルか）
        4. 文字数（最低3000文字以上か）
        5. セクションの完全性（すべてのセクションが含まれているか）
        6. 過去記事との連続性
        7. 伸びている記事の特徴を反映しているか
        
        レビュー結果はJSON形式で出力してください：
        {{
            "technical_accuracy": "技術的正確性の評価",
            "structure_quality": "構成の論理性の評価",
            "code_completeness": "コード例の完全性の評価",
            "word_count": 文字数,
            "section_completeness": "セクションの完全性の評価",
            "continuity": "過去記事との連続性の評価",
            "trending_features_reflection": "伸びている記事の特徴の反映度",
            "overall_score": 総合スコア（1-10）,
            "improvements": ["改善提案1", "改善提案2", ...],
            "approval": true/false
        }}
        """,
        agent=reviewer,
        expected_output="レビュー結果（JSON形式）。技術的正確性、構成、コード例、文字数、セクション完全性、連続性、改善提案、承認可否を含む。"
    )

