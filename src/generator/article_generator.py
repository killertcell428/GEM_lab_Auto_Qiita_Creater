"""記事生成モジュール - Google Gemini APIを使用して記事を生成"""
import google.generativeai as genai
from typing import Dict, Any, Optional, List
from src.config_loader import get_config
from src.prompt_loader import get_prompt_loader
from src.storage.qiita_items_manager import QiitaItemsManager


def initialize_gemini() -> None:
    """Gemini APIを初期化"""
    config = get_config()
    api_key = config.get("api", {}).get("gemini", {}).get("api_key", "")
    
    if not api_key:
        raise ValueError("GEMINI_API_KEYが設定されていません。.envファイルを確認してください。")
    
    genai.configure(api_key=api_key)


# 注意: この関数は古いワークフローで使用されていましたが、CrewAI版では使用されていません。
# 後方互換性のため残していますが、新しいコードでは使用しないでください。
def generate_article(
    research_data: Optional[Dict[str, Any]] = None,
    topic: Optional[str] = None,
    draft_content: Optional[str] = None,
    additional_requirements: Optional[str] = None
) -> Dict[str, Any]:
    """
    [非推奨] LLMを使用して記事を生成（ドラフト調整モード）
    
    この関数は古いワークフローで使用されていましたが、CrewAI版では使用されていません。
    新しいコードでは CrewAI の Orchestrator を使用してください。
    
    Args:
        research_data: リサーチデータ（後方互換性のため、使用されない）
        topic: トピック（後方互換性のため、使用されない）
        draft_content: ドラフト内容（Markdown形式）
        additional_requirements: 追加の要望や指示
        
    Returns:
        Dict[str, Any]: 生成された記事データ
            - title: str - 記事タイトル
            - content: str - 記事本文（Markdown形式）
            - tags: List[str] - タグのリスト
    """
    initialize_gemini()
    
    config = get_config()
    gemini_config = config.get("api", {}).get("gemini", {})
    model_name = gemini_config.get("model", "gemini-2.5-pro")
    # モデル名の正規化（models/プレフィックスを削除、genai.GenerativeModelはmodels/なしで受け取る）
    if model_name.startswith("models/"):
        model_name = model_name.replace("models/", "")
    temperature = gemini_config.get("temperature", 0.7)
    
    # ドラフト内容がない場合はエラー
    if not draft_content:
        raise ValueError("ドラフト内容が必要です。draft_contentパラメータを指定してください。")
    
    # プロンプトファイルから読み込む
    prompt_loader = get_prompt_loader()
    article_prompt_template = prompt_loader.load_prompt_or_default(
        "article_generation_prompt.md",
        """以下のドラフトをQiita用の記事に調整してください。

## ドラフト内容
{draft_content}

## 追加の要望
{additional_requirements}

要件:
- Qiitaの記事ルールに従う
- タイトルは明確で魅力的なものにする
- 技術的な内容を分かりやすく説明する
- コード例を含める（該当する場合）
- Markdown形式で記述する
- 見出し構造を適切に使用する
- 実践的で役立つ内容にする

出力形式:
タイトル: [記事のタイトル]

本文:
[Markdown形式の記事本文]

タグ: [関連する技術タグをカンマ区切りで3-5個]
"""
    )
    
    # 過去の投稿を取得
    qiita_items_manager = QiitaItemsManager()
    past_articles_summary = qiita_items_manager.get_items_summary(limit=10)
    
    # プロンプト改善機能を使用
    from src.generator.prompt_enhancer import PromptEnhancer
    enhancer = PromptEnhancer()
    
    # ベースプロンプトをフォーマット
    base_prompt = article_prompt_template.format(
        draft_content=draft_content,
        additional_requirements=additional_requirements or "（追加の要望なし）",
        past_articles_summary=past_articles_summary,
        trending_article_features="{trending_article_features}"  # プレースホルダー
    )
    
    # プロンプトを改善
    prompt = enhancer.enhance_article_prompt(base_prompt)
    
    try:
        # 高品質な長文生成のための設定（コストより品質重視）
        max_tokens = gemini_config.get("max_tokens", 200000)
        generation_config = genai.types.GenerationConfig(
            temperature=temperature,
            max_output_tokens=max_tokens,
        )
        
        model = genai.GenerativeModel(
            model_name=model_name,
            generation_config=generation_config
        )
        
        print(f"[INFO] 記事生成: モデル={model_name}, 最大トークン数={max_tokens}")
        
        response = model.generate_content(prompt)
        generated_text = response.text
        
        # レスポンスをパース
        article_data = parse_generated_article(generated_text)
        
        # 動的タグ生成（設定で有効な場合）
        if config.get("article", {}).get("dynamic_tags", True):
            dynamic_tags = generate_dynamic_tags(
                article_data.get("title", ""),
                article_data.get("content", "")
            )
            # デフォルトタグと動的タグをマージ
            default_tags = config.get("article", {}).get("default_tags", [])
            all_tags = list(set(default_tags + article_data.get("tags", []) + dynamic_tags))
            max_tags = config.get("article", {}).get("max_tags", 10)
            article_data["tags"] = all_tags[:max_tags]
        else:
            # 動的タグが無効な場合、デフォルトタグと生成されたタグをマージ
            default_tags = config.get("article", {}).get("default_tags", [])
            all_tags = list(set(default_tags + article_data.get("tags", [])))
            max_tags = config.get("article", {}).get("max_tags", 10)
            article_data["tags"] = all_tags[:max_tags]
        
        return article_data
    except Exception as e:
        raise Exception(f"記事の生成に失敗しました: {str(e)}")


def parse_generated_article(text: str) -> Dict[str, Any]:
    """
    生成されたテキストをパースして記事データに変換
    
    Args:
        text: LLMが生成したテキスト
        
    Returns:
        Dict[str, Any]: パースされた記事データ
    """
    lines = text.split("\n")
    title = ""
    content_lines = []
    tags = []
    
    # タイトル案セクションをスキップするフラグ
    skip_title_section = False
    found_main_title = False
    
    for i, line in enumerate(lines):
        # タイトル案セクションを検出してスキップ
        if "タイトル案" in line or (line.startswith("##") and "タイトル" in line):
            skip_title_section = True
            continue
        
        # タイトル案セクションの終了を検出
        if skip_title_section and (line.startswith("---") or (line.startswith("##") and "タイトル" not in line)):
            skip_title_section = False
        
        # タイトル案セクション中はスキップ
        if skip_title_section:
            continue
        
        # 最初の`# `で始まる行をタイトルとして取得（タイトル案セクションをスキップ後）
        if not found_main_title and line.startswith("# ") and "タイトル案" not in line:
            title = line.replace("# ", "").strip()
            found_main_title = True
            continue
        
        # タイトルが見つかった後の行は本文
        if found_main_title:
            # タグ行を検出
            if line.startswith("タグ:") or "タグ:" in line:
                tag_str = line.split("タグ:", 1)[1].strip() if "タグ:" in line else line.strip()
                tags = [tag.strip() for tag in tag_str.split(",") if tag.strip()]
            else:
                content_lines.append(line)
    
    # タイトルが取得できなかった場合、最初の`# `行を探す
    if not title:
        for line in lines:
            if line.startswith("# ") and "タイトル案" not in line:
                title = line.replace("# ", "").strip()
                break
    
    # デフォルト値の設定
    if not title:
        title = "生成された記事"
    
    content = "\n".join(content_lines).strip()
    if not content:
        content = text  # パースに失敗した場合は元のテキストを使用
    
    # デフォルトタグ
    config = get_config()
    default_tags = config.get("article", {}).get("default_tags", [])
    if not tags:
        tags = default_tags
    
    return {
        "title": title,
        "content": content,
        "tags": tags
    }


def generate_dynamic_tags(title: str, content: str) -> List[str]:
    """
    記事のタイトルと内容から動的にタグを生成
    
    Args:
        title: 記事のタイトル
        content: 記事の本文
        
    Returns:
        list[str]: 生成されたタグのリスト
    """
    config = get_config()
    gemini_config = config.get("api", {}).get("gemini", {})
    model_name = gemini_config.get("model", "gemini-2.5-pro")
    # モデル名の正規化（models/プレフィックスを削除、genai.GenerativeModelはmodels/なしで受け取る）
    if model_name.startswith("models/"):
        model_name = model_name.replace("models/", "")
    temperature = gemini_config.get("temperature", 0.7)
    
    # タグ生成用のプロンプト
    tag_prompt = f"""以下の記事のタイトルと内容を分析して、Qiitaに適したタグを生成してください。

タイトル: {title}

内容（抜粋）:
{content[:1000] if len(content) > 1000 else content}

## タグ生成の要件

以下のカテゴリから適切なタグを選んでください：

1. **言語・ツール系**: Python, R, Bash, Docker, Conda, JavaScript, TypeScript など
2. **目的・行為ベース**: 自動化, 環境構築, データ分析, 可視化, スクリプト, ワークフロー, CI/CD など
3. **プラットフォーム/インフラ系**: Linux, Ubuntu, Azure, GCP, AWS, GitHub Actions など
4. **ドメイン別**: バイオインフォマティクス, RNA-seq, 遺伝子発現解析, 機械学習, AI, データサイエンス など
5. **記事属性**: Tips, チュートリアル, 初心者向け, まとめ など

## 出力形式

以下の形式でタグを出力してください（カンマ区切り）：
タグ: タグ1, タグ2, タグ3, ...

重要:
- 記事の内容に最も関連するタグを優先する
- 6〜10個のタグを生成する
- Qiitaでよく使われるタグ名を使用する
- 重複を避ける
"""
    
    try:
        generation_config = genai.types.GenerationConfig(
            temperature=temperature,
            max_output_tokens=1000,  # タグ生成も詳細に
        )
        model = genai.GenerativeModel(
            model_name=model_name,
            generation_config=generation_config
        )
        
        response = model.generate_content(tag_prompt)
        response_text = response.text
        
        # タグをパース
        tags = []
        for line in response_text.split("\n"):
            if "タグ:" in line or "tags:" in line.lower():
                tag_str = line.split(":", 1)[1].strip() if ":" in line else line.strip()
                tags = [tag.strip() for tag in tag_str.split(",") if tag.strip()]
                break
        
        # タグが見つからない場合、最後の行を試す
        if not tags:
            last_line = response_text.strip().split("\n")[-1]
            if "," in last_line:
                tags = [tag.strip() for tag in last_line.split(",") if tag.strip()]
        
        # タグの正規化（空白削除、重複削除）
        tags = [tag.strip() for tag in tags if tag.strip()]
        tags = list(dict.fromkeys(tags))  # 重複を削除（順序を保持）
        
        return tags
    except Exception as e:
        # タグ生成に失敗した場合は空のリストを返す
        print(f"警告: 動的タグの生成に失敗しました: {str(e)}")
        return []

