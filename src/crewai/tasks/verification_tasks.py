"""Verification Tasks - Web検索・実装検証タスク"""
from crewai import Task
from typing import Dict, Any
from src.crewai.agents.web_researcher import create_web_researcher_agent
from src.crewai.agents.implementation_verifier import create_implementation_verifier_agent
from src.crewai.agents.knowledge_extractor import create_knowledge_extractor_agent
from src.crewai.state.article_state import ArticleState
import json


def create_web_research_task(topic: str, draft_content: str, context: Dict[str, Any]) -> Task:
    """
    Web検索タスクを作成
    
    Args:
        topic: 記事のトピック
        draft_content: ドラフト内容
        context: コンテキスト情報
    
    Returns:
        Web Research Task
    """
    web_researcher = create_web_researcher_agent()
    
    return Task(
        description=f"""
        以下のトピックについて、公式ドキュメントや実装例をWeb検索してください：
        
        トピック: {topic}
        
        ドラフト内容（抜粋）:
        {draft_content[:1000] if len(draft_content) > 1000 else draft_content}...
        
        検索項目:
        1. 公式ドキュメント（使用するライブラリの公式ドキュメント）
        2. 実装例ブログ（実際のコード例を含むブログ記事）
        3. GitHubリポジトリ（関連する実装例）
        4. ベストプラクティス（よく使われる手法やパターン）
        
        検索結果はJSON形式で出力してください：
        {{
            "official_docs": [
                {{
                    "title": "ドキュメントタイトル",
                    "url": "URL",
                    "snippet": "重要な情報の抜粋"
                }}
            ],
            "implementation_examples": [
                {{
                    "title": "記事タイトル",
                    "url": "URL",
                    "snippet": "コード例の抜粋"
                }}
            ],
            "best_practices": ["ベストプラクティス1", "ベストプラクティス2"],
            "summary": "検索結果の要約"
        }}
        """,
        agent=web_researcher,
        expected_output="Web検索結果（JSON形式）。公式ドキュメント、実装例、ベストプラクティスを含む。"
    )


def create_implementation_verification_task(
    code_snippets: str,
    web_research_results: str,
    context: Dict[str, Any]
) -> Task:
    """
    実装検証タスクを作成
    
    Args:
        code_snippets: 検証するコードスニペット
        web_research_results: Web検索結果
        context: コンテキスト情報
    
    Returns:
        Implementation Verification Task
    """
    verifier = create_implementation_verifier_agent()
    
    return Task(
        description=f"""
        以下のコードを実際に実行して動作確認し、実装ノウハウを抽出してください：
        
        検証するコード:
        {code_snippets}
        
        Web検索結果（参考）:
        {web_research_results}
        
        検証項目:
        1. コードが正常に実行できるか
        2. エラーが発生する場合、原因と修正案
        3. 実行時間やリソース使用量
        4. 実装過程でつまずいたポイント
        5. 最適化のヒント
        
        検証結果はJSON形式で出力してください：
        {{
            "verification_results": [
                {{
                    "code_snippet": "検証したコード",
                    "success": true/false,
                    "output": "実行結果",
                    "error": "エラーメッセージ（失敗時）",
                    "execution_time": 実行時間（秒）,
                    "notes": "実装時の注意点"
                }}
            ],
            "common_errors": [
                {{
                    "error_type": "エラーの種類",
                    "cause": "原因",
                    "solution": "解決方法"
                }}
            ],
            "optimization_tips": ["最適化のヒント1", "最適化のヒント2"],
            "implementation_notes": "実装過程で得たノウハウ"
        }}
        """,
        agent=verifier,
        expected_output="実装検証結果（JSON形式）。実行結果、エラー情報、最適化のヒント、実装ノウハウを含む。"
    )


def create_knowledge_extraction_task(
    verification_results: str,
    web_research_results: str,
    context: Dict[str, Any]
) -> Task:
    """
    ノウハウ抽出タスクを作成
    
    Args:
        verification_results: 実装検証結果
        web_research_results: Web検索結果
        context: コンテキスト情報
    
    Returns:
        Knowledge Extraction Task
    """
    extractor = create_knowledge_extractor_agent()
    
    return Task(
        description=f"""
        実装検証の過程で得られた知見から、読者に価値のあるノウハウを抽出してください：
        
        実装検証結果:
        {verification_results}
        
        Web検索結果:
        {web_research_results}
        
        抽出項目:
        1. 公式ドキュメントに書かれていない実践的な知見
        2. つまずきポイントと対処法
        3. パフォーマンス最適化のヒント
        4. よくあるエラーと解決方法
        5. 実務で使えるテクニック
        
        抽出結果はJSON形式で出力してください：
        {{
            "practical_insights": [
                {{
                    "insight": "実践的な知見",
                    "context": "どのような状況で役立つか",
                    "example": "具体例"
                }}
            ],
            "troubleshooting_guide": [
                {{
                    "problem": "問題の説明",
                    "solution": "解決方法",
                    "prevention": "予防策"
                }}
            ],
            "performance_tips": ["パフォーマンス最適化のヒント1", "ヒント2"],
            "best_practices": ["ベストプラクティス1", "ベストプラクティス2"],
            "summary": "ノウハウの要約"
        }}
        """,
        agent=extractor,
        expected_output="ノウハウ抽出結果（JSON形式）。実践的な知見、トラブルシューティングガイド、最適化のヒントを含む。"
    )

