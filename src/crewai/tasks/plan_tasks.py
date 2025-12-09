"""Plan Phase Tasks - 調査・構成設計タスク"""
from crewai import Task
from typing import Dict, Any
from src.crewai.agents.researcher import create_researcher_agent
from src.crewai.agents.planner import create_planner_agent
from src.crewai.agents.web_researcher import create_web_researcher_agent
from src.prompt_loader import get_prompt_loader


def create_research_task(topic: str, context: Dict[str, Any]) -> Task:
    """
    調査タスクを作成
    
    Args:
        topic: 調査するトピック
        context: コンテキスト情報（過去記事データなど）
    
    Returns:
        Research Task
    """
    researcher = create_researcher_agent()
    
    # プロンプトを読み込む
    prompt_loader = get_prompt_loader()
    research_prompt_template = prompt_loader.load_prompt_or_default(
        "research_prompt.md",
        """以下のトピックについて、技術記事のネタとして適切なリサーチを行ってください。

## トピック
{topic}

## 調査項目
1. 技術の概要・背景
2. 最新のトレンド・動向
3. 実装方法・ベストプラクティス
4. バイオインフォマティクス領域への応用可能性
5. 関連するツール・ライブラリ

調査結果はMarkdown形式でまとめてください。
"""
    )
    
    # 過去記事データを取得
    past_articles = context.get("past_articles", "なし")
    article_type_description = context.get("article_type_description", "")
    
    research_prompt = research_prompt_template.format(
        topic=topic,
        article_type_description=article_type_description
    )
    
    return Task(
        description=f"""
        以下のトピックについて深く調査してください：
        
        {research_prompt}
        
        過去の記事データ:
        {past_articles}
        
        調査結果は、後続のPlanner Agentが記事構成を設計する際に使用されます。
        実装可能な具体的な情報を含め、技術的な正確性を保ってください。
        """,
        agent=researcher,
        expected_output="調査レポート（Markdown形式）。技術の概要、最新トレンド、実装方法、バイオインフォマティクスへの応用可能性、関連ツールを含む。"
    )


def create_planning_task(research_report: str, context: Dict[str, Any]) -> Task:
    """
    構成設計タスクを作成
    
    Args:
        research_report: 調査レポート
        context: コンテキスト情報
    
    Returns:
        Planning Task
    """
    planner = create_planner_agent()
    
    past_articles = context.get("past_articles", "なし")
    article_type = context.get("article_type", "一般技術記事")
    trending_features = context.get("trending_features", "")
    
    return Task(
        description=f"""
        以下の調査結果を基に、Qiita記事の構成を設計してください：
        
        調査結果:
        {research_report}
        
        過去の記事:
        {past_articles}
        
        記事タイプ: {article_type}
        
        伸びている記事の特徴:
        {trending_features}
        
        設計項目:
        1. 記事タイトル（30-60文字、検索されやすいキーワードを含む）
        2. 対象読者（明確に定義）
        3. セクション構成（見出しと内容の概要）
           - 導入：背景・課題・なぜ重要か
           - トピックの概要
           - 技術的な仕組み・実装
           - 実践編：動くコード例
           - 応用・発展・実用例
           - まとめ
           - FAQ
           - 参考文献
        4. コード例の配置（どのセクションに何を配置するか）
        5. 過去記事との連続性（ストーリー性）
        6. 伸びている記事の特徴を参考にした改善点
        
        構成はJSON形式で出力してください。以下の構造で：
        {{
            "title": "記事タイトル",
            "target_audience": "対象読者",
            "sections": [
                {{
                    "heading": "見出し",
                    "content_outline": "内容の概要",
                    "code_examples": ["コード例の説明"]
                }}
            ],
            "code_placement": [
                {{
                    "section": "セクション名",
                    "description": "コード例の説明"
                }}
            ],
            "continuity": "過去記事との連続性の説明",
            "improvements": ["伸びている記事の特徴を参考にした改善点"]
        }}
        """,
        agent=planner,
        expected_output="記事構成プラン（JSON形式）。タイトル、対象読者、セクション構成、コード例の配置、過去記事との連続性、改善点を含む。"
    )

