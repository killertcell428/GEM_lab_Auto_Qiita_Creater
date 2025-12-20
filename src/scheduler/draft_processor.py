"""ドラフト記事の前処理（パス修正、画像・テーブル処理）"""
from pathlib import Path
from typing import Dict, Any, Optional, List
import re
import sys

# プロジェクトルートをパスに追加
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from src.publisher.github_image_uploader import GitHubImageUploader


class DraftProcessor:
    """ドラフト記事をQiita投稿用に前処理するクラス"""
    
    def __init__(self, draft_dir: Optional[Path] = None):
        """
        初期化
        
        Args:
            draft_dir: ドラフトファイルのディレクトリ（デフォルト: data/drafts/docs）
        """
        if draft_dir is None:
            self.draft_dir = project_root / "data" / "drafts" / "docs"
        else:
            self.draft_dir = Path(draft_dir)
        
        self.tables_dir = self.draft_dir / "blog" / "tables"
        self.plots_dir = self.draft_dir / "plots"
        self.slide_data_dir = self.draft_dir / "slide_data"
        self.image_uploader = GitHubImageUploader()
    
    def process_content(self, content: str, draft_file_path: Path) -> str:
        """
        記事内容を処理（パス修正、画像・テーブル処理）
        
        Args:
            content: 元の記事内容
            draft_file_path: ドラフトファイルのパス
            
        Returns:
            処理済みの記事内容
        """
        processed = content
        
        # 1. 相対パス参照を削除または修正
        processed = self._fix_relative_paths(processed, draft_file_path)
        
        # 2. テーブルファイルを記事内に埋め込む
        processed = self._embed_tables(processed)
        
        # 3. 画像参照を処理（外部URLに変換するか、参照を削除）
        processed = self._process_images(processed)
        
        # 4. コードブロック内のパス参照を処理
        processed = self._fix_code_block_paths(processed)
        
        return processed
    
    def _fix_relative_paths(self, content: str, draft_file_path: Path) -> str:
        """相対パス参照を修正"""
        # docs/ で始まる相対パス参照を削除（Qiitaでは参照できないため）
        # 例: [研究目的と背景](docs/research_objective.md) → [研究目的と背景](#)
        content = re.sub(
            r'\[([^\]]+)\]\(docs/[^\)]+\)',
            r'[\1](#)',
            content
        )
        
        # reference/ で始まるパス参照を削除
        content = re.sub(
            r'\[([^\]]+)\]\(reference/[^\)]+\)',
            r'[\1](#)',
            content
        )
        
        # simulation/ で始まるパス参照を削除
        content = re.sub(
            r'\[([^\]]+)\]\(simulation/[^\)]+\)',
            r'[\1](#)',
            content
        )
        
        # analysis/ で始まるパス参照を削除
        content = re.sub(
            r'\[([^\]]+)\]\(analysis/[^\)]+\)',
            r'[\1](#)',
            content
        )
        
        return content
    
    def _embed_tables(self, content: str) -> str:
        """テーブルファイルを記事内に埋め込む"""
        # ```csv:tables/XX_xxx.csv または ```csv:tables/XX_xxx.tsv の形式を検出
        pattern = r'```csv:tables/([^\s`]+)'
        
        def replace_table(match):
            table_file = match.group(1)
            table_path = self.tables_dir / table_file
            
            if not table_path.exists():
                print(f"[WARN] テーブルファイルが見つかりません: {table_path}")
                return match.group(0)  # 元のまま返す
            
            try:
                # CSVファイルを読み込んでMarkdownテーブル形式に変換
                with open(table_path, "r", encoding="utf-8") as f:
                    lines = f.readlines()
                
                if not lines:
                    return match.group(0)
                
                # ヘッダー行を取得
                header = lines[0].strip()
                if not header:
                    return match.group(0)
                
                # CSVをMarkdownテーブルに変換
                markdown_table = self._csv_to_markdown_table(lines)
                
                # 元のコードブロックを置き換え
                return f"```csv\n{markdown_table}\n```"
                
            except Exception as e:
                print(f"[WARN] テーブルファイルの読み込みエラー ({table_path}): {str(e)}")
                return match.group(0)
        
        content = re.sub(pattern, replace_table, content)
        return content
    
    def _csv_to_markdown_table(self, lines: List[str]) -> str:
        """CSV/TSVをMarkdownテーブル形式に変換"""
        if not lines:
            return ""
        
        import csv as csv_module
        from io import StringIO
        
        # ファイルの内容を結合
        content = "\n".join(lines)
        
        # 区切り文字を判定（TSVかCSVか）
        delimiter = "\t" if "\t" in content else ","
        
        # CSV/TSVをパース
        rows = []
        try:
            reader = csv_module.reader(StringIO(content), delimiter=delimiter)
            for row in reader:
                if row:  # 空行をスキップ
                    rows.append([cell.strip() for cell in row])
        except Exception as e:
            print(f"[WARN] CSV/TSVパースエラー: {str(e)}")
            # フォールバック: 簡易パース
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                if delimiter == "\t":
                    cells = [cell.strip() for cell in line.split("\t")]
                else:
                    cells = [cell.strip() for cell in line.split(",")]
                rows.append(cells)
        
        if not rows:
            return ""
        
        # 最大列数を取得
        max_cols = max(len(row) for row in rows) if rows else 0
        
        # すべての行を同じ列数に揃える
        for row in rows:
            while len(row) < max_cols:
                row.append("")
        
        # Markdownテーブル形式に変換
        markdown_lines = []
        
        # ヘッダー行
        if rows:
            header = rows[0]
            markdown_lines.append("| " + " | ".join(header) + " |")
            markdown_lines.append("| " + " | ".join(["---"] * len(header)) + " |")
            
            # データ行
            for row in rows[1:]:
                # セル内のパイプ文字をエスケープ
                escaped_row = [cell.replace("|", "\\|") for cell in row]
                markdown_lines.append("| " + " | ".join(escaped_row) + " |")
        
        return "\n".join(markdown_lines)
    
    def _process_images(self, content: str) -> str:
        """画像参照を処理（GitHubにアップロードしてURLを取得）"""
        # 画像ファイルへの参照を検出
        # 例: ![説明](volcano_sex_M_vs_F.png) または ![説明](analysis/out/plots/volcano_sex_M_vs_F.png)
        pattern = r'!\[([^\]]*)\]\(([^\)]+\.(png|jpg|jpeg|gif|svg|pdf))\)'
        
        def replace_image(match):
            alt_text = match.group(1)
            image_path = match.group(2)
            ext = match.group(3)
            
            # 画像ファイルのパスを解決
            image_file = Path(image_path).name  # ファイル名のみ取得
            
            # plotsディレクトリ内を検索
            possible_paths = [
                self.plots_dir / image_file,
                self.draft_dir / "blog" / image_file,
                self.draft_dir / image_file,
            ]
            
            image_file_path = None
            for path in possible_paths:
                if path.exists():
                    image_file_path = path
                    break
            
            if image_file_path:
                # 画像が見つかった場合、GitHubにアップロードを試みる
                print(f"[INFO] 画像ファイルが見つかりました: {image_file}")
                
                if self.image_uploader.is_configured():
                    # GitHubにアップロード
                    raw_url = self.image_uploader.upload_image(image_file_path)
                    if raw_url:
                        # アップロード成功：画像参照を更新
                        print(f"[INFO] 画像をGitHubにアップロードしました: {raw_url}")
                        return f"![{alt_text}]({raw_url})"
                    else:
                        # アップロード失敗：コメントアウト
                        print(f"[WARN] 画像のアップロードに失敗しました: {image_file}")
                        return f"<!-- TODO: 画像URLを設定してください: {image_file} -->\n<!-- ![{alt_text}](画像URL) -->"
                else:
                    # GitHub設定がない場合：コメントアウト
                    print(f"[INFO] GitHub設定がないため、画像URLを手動で設定してください: {image_file}")
                    return f"<!-- TODO: 画像URLを設定してください: {image_file} -->\n<!-- ![{alt_text}](画像URL) -->"
            else:
                print(f"[WARN] 画像ファイルが見つかりません: {image_path}")
                # 画像参照をコメントアウト
                return f"<!-- 画像参照を削除: {image_path} -->"
        
        content = re.sub(pattern, replace_image, content)
        return content
    
    def _fix_code_block_paths(self, content: str) -> str:
        """コードブロック内のパス参照を修正"""
        # コードブロック内のパスは基本的にそのまま残す
        # ただし、明らかにローカルパスのみの参照はコメントアウト
        
        # コードブロック内のパス参照を検出して処理
        # 例: ```bash:code_snippets/01_backbone_download.sh
        pattern = r'```(\w+):([^\s`]+)'
        
        def replace_code_path(match):
            lang = match.group(1)
            code_path = match.group(2)
            
            # コードブロックのパス指定を削除（Qiitaでは使えない）
            return f"```{lang}"
        
        content = re.sub(pattern, replace_code_path, content)
        return content
