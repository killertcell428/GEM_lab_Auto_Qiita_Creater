"""プロンプト改善モジュール - 伸びている記事の特徴をプロンプトに反映"""
from typing import Dict, Any, Optional
from src.analyzer.article_analyzer import ArticleAnalyzer
from src.analyzer.feedback_collector import FeedbackCollector


class PromptEnhancer:
    """プロンプトを改善するクラス"""
    
    def __init__(self):
        """初期化"""
        self.analyzer = ArticleAnalyzer()
        self.feedback_collector = FeedbackCollector()
    
    def get_enhanced_prompt_context(self) -> str:
        """
        改善されたプロンプトコンテキストを取得
        
        Returns:
            str: プロンプトに追加するコンテキスト
        """
        # 保存された記事特徴を読み込む
        features = self.analyzer.load_features()
        
        if not features:
            return self._get_default_context()
        
        # 高パフォーマンス記事の情報を取得
        high_performance = self.feedback_collector.get_high_performance_articles(
            min_likes=3,
            min_views=50
        )
        
        context_parts = []
        
        # 統計情報
        context_parts.append("## 伸びている記事の統計情報")
        context_parts.append(f"- 平均タイトル文字数: {features.get('avg_title_length', 0)}文字")
        context_parts.append(f"- 平均本文文字数: {features.get('avg_content_length', 0)}文字")
        context_parts.append(f"- 平均コードブロック数: {features.get('avg_code_blocks', 0)}個")
        context_parts.append(f"- 平均見出し数: {features.get('avg_headings', 0)}個")
        context_parts.append(f"- 平均タグ数: {features.get('avg_tags', 0)}個")
        context_parts.append(f"- コード例を含む記事の割合: {features.get('has_code_ratio', 0) * 100:.0f}%")
        
        # 共通の構造
        if features.get("common_structure"):
            context_parts.append("")
            context_parts.append("## 効果的な記事構造")
            context_parts.append(features.get("common_structure"))
        
        # 重要なポイント
        if features.get("common_key_points"):
            context_parts.append("")
            context_parts.append("## 読者が価値を感じるポイント")
            for i, point in enumerate(features.get("common_key_points", []), 1):
                context_parts.append(f"{i}. {point}")
        
        # 高パフォーマンス記事の例
        if high_performance:
            context_parts.append("")
            context_parts.append("## 高パフォーマンス記事の例")
            for i, article in enumerate(high_performance[:3], 1):
                title = article.get("title", "")
                likes = article.get("likes_count", 0)
                views = article.get("page_views_count", 0)
                context_parts.append(f"{i}. {title} (いいね: {likes}, 閲覧: {views})")
        
        return "\n".join(context_parts)
    
    def _get_default_context(self) -> str:
        """
        デフォルトのコンテキストを取得（分析データがない場合）
        
        Returns:
            str: デフォルトコンテキスト
        """
        return """## 伸びている記事の一般的な特徴

- タイトルは30-60文字程度で、具体的で検索されやすいキーワードを含む
- 本文は3000文字以上で、実践的な内容を含む
- コード例を含む（特に実装例やサンプルコード）
- 見出し構造が明確で、読みやすい構成
- タグは5-10個程度で、適切にカテゴライズされている
- 読者が実際に試せる内容
- 技術的な正確性と実践性のバランスが取れている"""
    
    def enhance_article_prompt(self, base_prompt: str) -> str:
        """
        記事生成プロンプトを改善
        
        Args:
            base_prompt: ベースとなるプロンプト
            
        Returns:
            str: 改善されたプロンプト
        """
        enhanced_context = self.get_enhanced_prompt_context()
        
        # プロンプトに改善コンテキストを追加
        if "{trending_article_features}" in base_prompt:
            enhanced_prompt = base_prompt.replace("{trending_article_features}", enhanced_context)
        else:
            # プロンプトの最後に追加
            enhanced_prompt = f"""{base_prompt}

---

{enhanced_context}
"""
        
        return enhanced_prompt

