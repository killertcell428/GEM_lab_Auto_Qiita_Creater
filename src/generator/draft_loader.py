"""ドラフト読み込みモジュール - ファイルからドラフトを読み込む"""
from pathlib import Path
from typing import Dict, Any, Optional
import re


def load_draft_from_file(file_path: str) -> Dict[str, Any]:
    """
    ファイルからドラフトを読み込む
    
    Args:
        file_path: ドラフトファイルのパス（Markdown形式）
        
    Returns:
        Dict[str, Any]: ドラフトデータ
            - title: Optional[str] - タイトル（見つかった場合）
            - content: str - 本文
            - tags: List[str] - タグ（見つかった場合）
            - raw_content: str - 生のファイル内容
    """
    draft_path = Path(file_path)
    
    if not draft_path.exists():
        raise FileNotFoundError(f"ドラフトファイルが見つかりません: {file_path}")
    
    if not draft_path.is_file():
        raise ValueError(f"指定されたパスはファイルではありません: {file_path}")
    
    # ファイルを読み込む
    try:
        with open(draft_path, "r", encoding="utf-8") as f:
            raw_content = f.read()
    except Exception as e:
        raise Exception(f"ドラフトファイルの読み込みに失敗しました: {str(e)}")
    
    # ドラフトをパース
    parsed = parse_draft(raw_content)
    parsed["raw_content"] = raw_content
    parsed["file_path"] = str(draft_path)
    
    return parsed


def parse_draft(content: str) -> Dict[str, Any]:
    """
    Markdown形式のドラフトをパース
    
    Args:
        content: Markdown形式のドラフト内容
        
    Returns:
        Dict[str, Any]: パースされたドラフトデータ
            - title: Optional[str] - タイトル
            - content: str - 本文
            - tags: List[str] - タグ
    """
    lines = content.split("\n")
    title = None
    tags = []
    content_lines = []
    
    # 最初の `# ` で始まる行をタイトルとして取得
    for i, line in enumerate(lines):
        if line.startswith("# "):
            title = line.replace("# ", "").strip()
            # タイトル以降の行を本文として扱う
            content_lines = lines[i + 1:]
            break
        elif line.strip():
            content_lines.append(line)
    
    # タグを検出（`タグ:` または `tags:` で始まる行）
    for i, line in enumerate(content_lines):
        if re.match(r"^タグ[:：]\s*", line, re.IGNORECASE) or re.match(r"^tags[:：]\s*", line, re.IGNORECASE):
            tag_str = re.split(r"[:：]\s*", line, 1)[1] if ":" in line or "：" in line else line.strip()
            tags = [tag.strip() for tag in re.split(r"[,，、]", tag_str) if tag.strip()]
            # タグ行を本文から削除
            content_lines.pop(i)
            break
    
    # 本文を結合
    content = "\n".join(content_lines).strip()
    
    # タイトルが見つからなかった場合、最初の行をタイトル候補とする
    if not title and content_lines:
        first_line = content_lines[0].strip()
        if len(first_line) < 100 and not first_line.startswith("#"):
            title = first_line
            content = "\n".join(content_lines[1:]).strip()
    
    return {
        "title": title,
        "content": content,
        "tags": tags
    }

