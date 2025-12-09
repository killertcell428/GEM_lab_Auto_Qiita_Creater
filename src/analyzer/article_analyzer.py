"""記事特徴分析モジュール - 伸びている記事の特徴を分析"""
import google.generativeai as genai
from typing import Dict, Any, List, Optional
from pathlib import Path
import json
from datetime import datetime
from src.config_loader import get_config


class ArticleAnalyzer:
    """記事の特徴を分析するクラス"""
    
    def __init__(self):
        """初期化"""
        config = get_config()
        gemini_config = config.get("api", {}).get("gemini", {})
        self.api_key = gemini_config.get("api_key", "")
        
        if self.api_key:
            genai.configure(api_key=self.api_key)
            model_name = gemini_config.get("model", "gemini-2.5-pro")
            if model_name.startswith("models/"):
                model_name = model_name.replace("models/", "")
            self.model = genai.GenerativeModel(model_name)
        else:
            self.model = None
        
        # データ保存ディレクトリ
        self.data_dir = Path("data/analysis")
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.features_file = self.data_dir / "article_features.json"
    
    def analyze_article(self, article: Dict[str, Any]) -> Dict[str, Any]:
        """
        単一の記事を分析
        
        Args:
            article: 記事データ（Qiita APIのレスポンス形式）
            
        Returns:
            Dict[str, Any]: 分析結果
                - title_length: int - タイトル文字数
                - content_length: int - 本文文字数
                - has_code: bool - コード例の有無
                - code_blocks_count: int - コードブロックの数
                - heading_count: int - 見出しの数
                - tags_count: int - タグ数
                - likes_count: int - いいね数
                - page_views_count: int - 閲覧数
                - structure: str - 記事の構造（LLM分析結果）
                - key_points: List[str] - 重要なポイント（LLM抽出）
        """
        title = article.get("title", "")
        body = article.get("body", "")
        tags = article.get("tags", [])
        likes_count = article.get("likes_count", 0)
        page_views_count = article.get("page_views_count", 0)
        
        # 基本的な統計
        features = {
            "title_length": len(title),
            "content_length": len(body),
            "has_code": "```" in body or "`" in body,
            "code_blocks_count": body.count("```"),
            "heading_count": body.count("#"),
            "tags_count": len(tags),
            "likes_count": likes_count,
            "page_views_count": page_views_count,
            "article_id": article.get("id", ""),
            "url": article.get("url", ""),
            "created_at": article.get("created_at", ""),
        }
        
        # LLMで詳細分析（APIキーがある場合）
        if self.model and body:
            try:
                llm_analysis = self._analyze_with_llm(title, body, likes_count, page_views_count)
                features.update(llm_analysis)
            except Exception as e:
                print(f"[WARN] LLM分析に失敗: {str(e)}")
        
        return features
    
    def _analyze_with_llm(self, title: str, body: str, likes_count: int, page_views_count: int) -> Dict[str, Any]:
        """
        LLMを使って記事の特徴を分析
        
        Args:
            title: 記事タイトル
            body: 記事本文
            likes_count: いいね数
            page_views_count: 閲覧数
            
        Returns:
            Dict[str, Any]: LLM分析結果
        """
        prompt = f"""以下のQiita記事を分析して、なぜこの記事が人気（いいね数: {likes_count}, 閲覧数: {page_views_count}）なのか、その特徴を抽出してください。

## タイトル
{title}

## 本文（抜粋）
{body[:2000] if len(body) > 2000 else body}

## 分析項目
1. **記事の構造**: どのような構成になっているか（導入→説明→実装→まとめなど）
2. **重要なポイント**: 読者が価値を感じるポイントを3-5個抽出
3. **技術的な深さ**: 初心者向けか、上級者向けか
4. **実践性**: 実際に使える内容かどうか
5. **タイトルの特徴**: なぜこのタイトルが効果的か

## 出力形式
以下のJSON形式で出力してください：
{{
  "structure": "記事の構造の説明",
  "key_points": ["ポイント1", "ポイント2", "ポイント3"],
  "technical_depth": "初心者向け/中級者向け/上級者向け",
  "practicality": "実践性の説明",
  "title_effectiveness": "タイトルの効果的な理由"
}}
"""
        
        try:
            response = self.model.generate_content(prompt)
            response_text = response.text.strip()
            
            # JSONを抽出（コードブロックがある場合）
            if "```json" in response_text:
                response_text = response_text.split("```json")[1].split("```")[0].strip()
            elif "```" in response_text:
                response_text = response_text.split("```")[1].split("```")[0].strip()
            
            # JSONをパース
            analysis = json.loads(response_text)
            return analysis
            
        except Exception as e:
            print(f"[WARN] LLM分析のパースに失敗: {str(e)}")
            return {}
    
    def analyze_multiple_articles(self, articles: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        複数の記事を分析
        
        Args:
            articles: 記事データのリスト
            
        Returns:
            List[Dict[str, Any]]: 分析結果のリスト
        """
        analyzed = []
        for i, article in enumerate(articles, 1):
            print(f"[INFO] 記事を分析中 ({i}/{len(articles)})...")
            features = self.analyze_article(article)
            analyzed.append(features)
        
        return analyzed
    
    def extract_common_features(self, analyzed_articles: List[Dict[str, Any]], top_n: int = 10) -> Dict[str, Any]:
        """
        分析結果から共通の特徴を抽出
        
        Args:
            analyzed_articles: 分析結果のリスト
            top_n: 上位N件の記事を対象にする
            
        Returns:
            Dict[str, Any]: 共通特徴
        """
        # いいね数と閲覧数の合計でソート
        sorted_articles = sorted(
            analyzed_articles,
            key=lambda x: x.get("likes_count", 0) + (x.get("page_views_count") or 0) * 0.1,
            reverse=True
        )
        
        top_articles = sorted_articles[:top_n]
        
        if not top_articles:
            return {}
        
        # 統計を計算
        avg_title_length = sum(a.get("title_length", 0) for a in top_articles) / len(top_articles)
        avg_content_length = sum(a.get("content_length", 0) for a in top_articles) / len(top_articles)
        avg_code_blocks = sum(a.get("code_blocks_count", 0) for a in top_articles) / len(top_articles)
        avg_headings = sum(a.get("heading_count", 0) for a in top_articles) / len(top_articles)
        avg_tags = sum(a.get("tags_count", 0) for a in top_articles) / len(top_articles)
        has_code_ratio = sum(1 for a in top_articles if a.get("has_code", False)) / len(top_articles)
        
        # 構造とポイントを集約
        structures = [a.get("structure", "") for a in top_articles if a.get("structure")]
        key_points = []
        for a in top_articles:
            points = a.get("key_points", [])
            if isinstance(points, list):
                key_points.extend(points)
        
        # 最も多い構造とポイントを抽出
        from collections import Counter
        structure_counter = Counter(structures)
        most_common_structure = structure_counter.most_common(1)[0][0] if structure_counter else ""
        
        point_counter = Counter(key_points)
        most_common_points = [point for point, _ in point_counter.most_common(5)]
        
        return {
            "avg_title_length": int(avg_title_length),
            "avg_content_length": int(avg_content_length),
            "avg_code_blocks": round(avg_code_blocks, 1),
            "avg_headings": int(avg_headings),
            "avg_tags": round(avg_tags, 1),
            "has_code_ratio": round(has_code_ratio, 2),
            "common_structure": most_common_structure,
            "common_key_points": most_common_points,
            "analyzed_count": len(top_articles),
            "analyzed_at": datetime.now().isoformat()
        }
    
    def save_features(self, features: Dict[str, Any]) -> None:
        """
        分析結果をJSONファイルに保存
        
        Args:
            features: 保存する特徴データ
        """
        data = {
            "updated_at": datetime.now().isoformat(),
            "features": features
        }
        
        try:
            with open(self.features_file, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, indent=2)
            print(f"[OK] 記事特徴を保存しました: {self.features_file}")
        except Exception as e:
            print(f"[WARN] 記事特徴の保存に失敗: {str(e)}")
    
    def load_features(self) -> Optional[Dict[str, Any]]:
        """
        保存された記事特徴を読み込む
        
        Returns:
            Optional[Dict[str, Any]]: 記事特徴データ（存在しない場合はNone）
        """
        if not self.features_file.exists():
            return None
        
        try:
            with open(self.features_file, "r", encoding="utf-8") as f:
                data = json.load(f)
            return data.get("features")
        except Exception as e:
            print(f"[WARN] 記事特徴の読み込みに失敗: {str(e)}")
            return None

