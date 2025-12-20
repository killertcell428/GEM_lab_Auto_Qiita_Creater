# 第3回：データ準備「バックボーンゲノムの取得と染色体リネーム」

## はじめに

前回は、ターミナル操作とconda環境の構築について解説しました。今回は、**バックボーンゲノムの取得とデータ準備**について詳しく説明します。

RNA-seq解析には、リファレンスゲノムが必要です。このプロジェクトでは、ニワトリ（Gallus gallus）のゲノムをバックボーンとして使用します。

---

## パイプライン全体図

この記事は、以下のパイプラインの**ステップ3：データ準備**に対応しています：

```
[1] イントロダクション
[2] 環境構築
[3] データ準備 ← 現在ここ
[4] シミュレーション
[5] アライメント・カウント
[6] 発現解析
[7] 可視化
[8] CNV解析
[9] まとめ
```

<!-- TODO: 画像URLを設定してください: パイプライン全体図（ステップ3をハイライト） -->
<!-- ![パイプライン全体図](画像URL) -->

---

## 前回の振り返り

**第2回で学んだこと**：
- ターミナルの基本操作（`pwd`, `ls`, `cd`, `mkdir`など）
- conda環境の作成と管理（`conda create`, `conda activate`）
- ディレクトリ構成の設計とシンボリックリンクの活用

**今回から始めること**：
- Ensemblからゲノムデータをダウンロード
- FASTA/GTFファイルの理解
- 染色体リネームの実行

---

**この記事で学べること**：
- Ensemblからゲノムデータをダウンロードする方法
- FASTA/GTFファイルとは何か
- 染色体リネームの方法
- ダウンロードしたデータの確認方法

**前提知識**：
- 前回の記事（環境構築）を完了していること
- ターミナルの基本的な操作ができること

---

## バックボーンゲノムの取得

### なぜニワトリ（Gallus gallus）なのか？

リオス科は「飛竜種」という設定で、現実世界では：

- **鳥類**：飛翔能力、性染色体（ZW型）、猛禽類の性差
- **爬虫類**：体軸伸長（Hox遺伝子のシフト）、尾の形態多様性

ニワトリゲノムを選んだ理由：

1. **鳥類の代表種**：飛翔に関わる遺伝子セットが揃っている
2. **ゲノム品質が高い**：Ensembl/NCBIで整備されたアノテーション（GRCg7b）
3. **性染色体が明確**：Z/W染色体があり、性決定様式（ZW型）がリオス科設定と一致
4. **比較ゲノム学でよく使われる**：他の鳥類・爬虫類との比較データが豊富

**参考記事**：
- [ゲノムデータベースの使い方（Ensembl/NCBI）](https://qiita.com/tags/genome) - Qiitaのゲノムデータベース記事
- [Ensembl公式サイト](https://www.ensembl.org/) - ゲノムデータベース
- [NCBI公式サイト](https://www.ncbi.nlm.nih.gov/) - ゲノムデータベース

<!-- TODO: 画像URLを設定してください: ニワトリゲノム選択の理由図 -->
<!-- ![ニワトリゲノム選択の理由](画像URL) -->

### Ensemblとは？

**Ensembl**は、様々な生物種のゲノム配列とアノテーションを提供するデータベースです。

**Ensemblの特徴**：
- 100以上の生物種のゲノムデータを提供
- 定期的に更新される（年2-3回）
- ゲノム配列（FASTA）とアノテーション（GTF/GFF）を提供

**FASTAファイルとは？**
- DNA/RNA/タンパク質配列を記録する標準的な形式
- ヘッダー行（`>`で始まる）と配列行で構成

**GTFファイルとは？**
- 遺伝子アノテーション（遺伝子の位置、エクソン、イントロンなど）を記録する形式
- 9列のタブ区切り形式

**参考記事**：
- [FASTA形式の説明](https://qiita.com/tags/fasta) - QiitaのFASTA記事
- [GTF/GFF形式の説明](https://qiita.com/tags/gtf) - QiitaのGTF記事
- [NCBI FASTA形式の説明](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) - NCBI公式

### Ensemblからのダウンロード

**ダウンロードの準備**：

```bash
# ダウンロード先のディレクトリに移動
cd /home/uecha/data/reference/rios

# ディレクトリが存在しない場合は作成
mkdir -p /home/uecha/data/reference/rios
```

**手動ダウンロード（初心者向け）**：

1. [Ensembl FTPサイト](https://ftp.ensembl.org/pub/)にアクセス
2. `release-112/fasta/gallus_gallus/dna/`に移動
3. `Gallus_gallus.GRCg7b.dna.primary_assembly.fa.gz`をダウンロード
4. 同様に、`release-112/gtf/gallus_gallus/`からGTFファイルもダウンロード

**手動ダウンロードのメリット**：
- ブラウザで確認しながらダウンロードできる
- ファイル名を確認してからダウンロードできる

**手動ダウンロードのデメリット**：
- 時間がかかる
- ファイル名が長くて間違えやすい

**スクリプトでのダウンロード（推奨）**：

```bash
# ダウンロードスクリプトを作成
cat > download_backbone.sh << 'EOF'
#!/usr/bin/env bash
set -euo pipefail

ENSEMBL_RELEASE="112"
SPECIES="gallus_gallus"
BASE_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}"

# 複数の候補URLを用意（ファイル名の違いに対応）
FA_CANDIDATES=(
  "${BASE_URL}/fasta/${SPECIES}/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz"
  "${BASE_URL}/fasta/${SPECIES}/dna/Gallus_gallus.GRCg7b.dna.toplevel.fa.gz"
  "${BASE_URL}/fasta/${SPECIES}/dna/Gallus_gallus.GRCg7b.dna.primary_assembly.fa.gz"
)

GTF_CANDIDATES=(
  "${BASE_URL}/gtf/${SPECIES}/Gallus_gallus.GRCg7b.112.gtf.gz"
  "${BASE_URL}/gtf/${SPECIES}/Gallus_gallus.GRCg7b.${ENSEMBL_RELEASE}.gtf.gz"
)

# 最初に利用可能なURLからダウンロード
download_first_available() {
  local outfile="$1"; shift
  for url in "$@"; do
    echo "[info] try ${url}"
    if wget -c "$url" -O "$outfile"; then
      return 0
    fi
  done
  echo "[error] all download attempts failed" >&2
  return 1
}

# FASTAファイルのダウンロード
download_first_available rios_ancestor.fa.gz "${FA_CANDIDATES[@]}"

# GTFファイルのダウンロード
download_first_available rios_ancestor.gtf.gz "${GTF_CANDIDATES[@]}"

# 解凍
gzip -d -f rios_ancestor.fa.gz
gzip -d -f rios_ancestor.gtf.gz

echo "[done] Download completed!"
EOF

# 実行権限を付与
chmod +x download_backbone.sh

# スクリプトを実行
bash download_backbone.sh
```

**スクリプトの説明**：
- `#!/usr/bin/env bash`：bashスクリプトであることを示す
- `set -euo pipefail`：エラー時に即座に停止
- `wget -c`：中断された場合に再開可能（`-c`：継続ダウンロード）
- `gzip -d`：gzipファイルを解凍

**ダウンロード時間**：
- ゲノムFASTAファイル：約1GB、数分〜数十分（ネットワーク速度による）
- GTFファイル：約100MB、数分

**ダウンロードが途中で止まった場合**：
```bash
# wgetの再開機能を使う（-cオプション）
wget -c https://...

# または、curlを使う
curl -C - -O https://...
```

**参考記事**：
- [wget/curlでファイルをダウンロードする方法](https://qiita.com/tags/wget) - Qiitaのダウンロード記事
- [bashスクリプトの書き方（初心者向け）](https://qiita.com/tags/bash) - Qiitaのbash記事

### ダウンロードしたファイルの確認

```bash
# ファイルサイズの確認
ls -lh /home/uecha/data/reference/rios/
# 出力例:
# -rw-r--r-- 1 uecha uecha 1.0G Dec 11 19:25 rios_ancestor.fa
# -rw-r--r-- 1 uecha uecha 100M Dec 11 19:30 rios_ancestor.gtf

# FASTAファイルの最初の数行を確認
head -5 /home/uecha/data/reference/rios/rios_ancestor.fa
# 出力例:
# >chr1 dna:primary_assembly ...
# ATGCGATCGATCG...
# ATGCGATCGATCG...

# GTFファイルの最初の数行を確認
head -5 /home/uecha/data/reference/rios/rios_ancestor.gtf
# 出力例:
# #!genome-build Gallus_gallus-7.0
# #!genome-version GRCg7b
# chr1	ensembl	gene	...
```

**コマンドの説明**：
- `ls -lh`：ファイルの詳細情報を表示（`-l`：詳細、`-h`：人間が読みやすい形式）
- `head -5`：ファイルの最初の5行を表示

**参考記事**：
- [Linuxコマンド一覧](https://eng-entrance.com/linux-command) - コマンドの解説サイト

---

## 染色体リネーム

### なぜリネームするのか？

`chr1` → `Rios1` などにリネームしたのは：

- **「架空種のゲノム」として扱いやすくする**ため
- 解析ツール（STAR, HISAT2, featureCounts）は染色体名を識別するため、一貫性が必要
- 将来的に「リオス科 vs ニワトリ」の比較解析をする際、混同を避けるため

**リネームの例**：
- `chr1` → `Rios1`
- `chr2` → `Rios2`
- `chrZ` → `RiosZ`
- `chrW` → `RiosW`
- `chrMT` → `RiosMT`

<!-- TODO: 画像URLを設定してください: 染色体リネームの概念図 -->
<!-- ![染色体リネームの概念図](画像URL) -->

### リネーム対応表の作成

まず、リネーム対応表（TSV形式）を作成します：

```bash
# リネーム対応表を作成
cat > chrom_rename.tsv << 'EOF'
old	new
chr1	Rios1
chr2	Rios2
chr3	Rios3
chr4	Rios4
chr5	Rios5
chr6	Rios6
chr7	Rios7
chr8	Rios8
chr9	Rios9
chr10	Rios10
chr11	Rios11
chr12	Rios12
chr13	Rios13
chr14	Rios14
chr15	Rios15
chr16	Rios16
chr17	Rios17
chr18	Rios18
chr19	Rios19
chr20	Rios20
chr21	Rios21
chr22	Rios22
chr23	Rios23
chr24	Rios24
chr25	Rios25
chr26	Rios26
chr27	Rios27
chr28	Rios28
chrZ	RiosZ
chrW	RiosW
chrMT	RiosMT
EOF
```

**TSV形式とは？**
- Tab Separated Values（タブ区切り値）の略
- 表形式のデータを保存する形式
- CSV（カンマ区切り）と似ているが、タブで区切る

### リネームスクリプトの実装

```python
#!/usr/bin/env python3
"""染色体名をリネームするスクリプト。

Usage:
  python scripts/rename_chrom.py input.fa chrom_rename.tsv output.fa
  python scripts/rename_chrom.py input.gtf chrom_rename.tsv output.gtf
"""

import sys
from pathlib import Path

def load_map(tsv_path: Path):
    """リネーム対応表を読み込む"""
    mapping = {}
    with tsv_path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("old"):
                continue
            old, new = line.rstrip().split("\t")
            mapping[old] = new
    return mapping

def rename_fasta(inp: Path, out: Path, mapping: dict):
    """FASTAファイルのヘッダーをリネーム"""
    with inp.open() as fin, out.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                rest = line[len(name)+1:]
                fout.write(f">{mapping.get(name, name)}{rest}")
            else:
                fout.write(line)

def rename_gtf(inp: Path, out: Path, mapping: dict):
    """GTFファイルの1列目（染色体名）をリネーム"""
    with inp.open() as fin, out.open("w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                fout.write(line)
                continue
            old = parts[0]
            parts[0] = mapping.get(old, old)
            fout.write("\t".join(parts) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 rename_chrom.py input_file rename_table.tsv output_file")
        sys.exit(1)
    
    inp_path = Path(sys.argv[1])
    map_path = Path(sys.argv[2])
    out_path = Path(sys.argv[3])
    
    mapping = load_map(map_path)
    
    if inp_path.suffix in [".fa", ".fasta"]:
        rename_fasta(inp_path, out_path, mapping)
    elif inp_path.suffix == ".gtf":
        rename_gtf(inp_path, out_path, mapping)
    else:
        print(f"Error: Unknown file type: {inp_path.suffix}")
        sys.exit(1)
    
    print(f"Renamed: {inp_path} -> {out_path}")
```

**スクリプトの説明**：
- `load_map`：TSVファイルからリネーム対応表を読み込む
- `rename_fasta`：FASTAファイルのヘッダー（`>`で始まる行）をリネーム
- `rename_gtf`：GTFファイルの1列目（染色体名）をリネーム

**Pythonスクリプトの実行方法**：

```bash
# Pythonがインストールされているか確認
python3 --version
# 出力例: Python 3.9.7

# スクリプトを保存（scripts/rename_chrom.py）
mkdir -p scripts
# 上記のPythonコードを scripts/rename_chrom.py に保存

# 実行権限を付与
chmod +x scripts/rename_chrom.py

# FASTAファイルのリネーム
python3 scripts/rename_chrom.py \
  /home/uecha/data/reference/rios/rios_ancestor.fa \
  chrom_rename.tsv \
  /home/uecha/data/reference/rios/rios_genome.fa

# GTFファイルのリネーム
python3 scripts/rename_chrom.py \
  /home/uecha/data/reference/rios/rios_ancestor.gtf \
  chrom_rename.tsv \
  /home/uecha/data/reference/rios/rios.gtf
```

**参考記事**：
- [Python入門（初心者向け）](https://qiita.com/tags/python) - QiitaのPython記事
- [Python公式チュートリアル](https://docs.python.org/ja/3/tutorial/) - 公式ドキュメント
- [ファイル操作の基本（Python）](https://qiita.com/tags/python) - Qiitaのファイル操作記事

### リネーム結果の確認

```bash
# リネームが正しく行われたか確認
head -5 /home/uecha/data/reference/rios/rios_genome.fa
# 出力例:
# >Rios1 dna:primary_assembly ...
# >Rios2 dna:primary_assembly ...

# 染色体数の確認
grep "^>Rios" /home/uecha/data/reference/rios/rios_genome.fa | wc -l
# 出力例: 50

# GTFファイルの確認
head -5 /home/uecha/data/reference/rios/rios.gtf
# 出力例:
# #!genome-build Gallus_gallus-7.0
# #!genome-version GRCg7b
# Rios1	ensembl	gene	...
```

**コマンドの説明**：
- `grep "^>Rios"`：`>`で始まり`Rios`を含む行を検索（FASTAヘッダーを検索）
- `wc -l`：行数をカウント

**参考記事**：
- [grepコマンドの使い方](https://qiita.com/tags/grep) - Qiitaのgrep記事

---

## シンボリックリンクの作成

前回作成したディレクトリ構成に、ゲノムファイルへのシンボリックリンクを作成します。

```bash
# プロジェクトディレクトリに移動
cd /home/uecha/work/MonsterGenome_Rios

# ゲノムファイルへのリンク作成
ln -sf /home/uecha/data/reference/rios/rios_genome.fa \
       /home/uecha/work/MonsterGenome_Rios/reference/rios_genome.fa

# GTFファイルへのリンク作成
ln -sf /home/uecha/data/reference/rios/rios.gtf \
       /home/uecha/work/MonsterGenome_Rios/reference/rios.gtf
```

**確認方法**：
```bash
# リンクが正しく作成されたか確認
ls -lh /home/uecha/work/MonsterGenome_Rios/reference/
# 出力例:
# lrwxrwxrwx 1 uecha uecha 45 Dec 11 19:28 rios_genome.fa -> /home/uecha/data/reference/rios/rios_genome.fa
# lrwxrwxrwx 1 uecha uecha 42 Dec 11 19:28 rios.gtf -> /home/uecha/data/reference/rios/rios.gtf
```

**出力の見方**：
- `lrwxrwxrwx`：最初の`l`はシンボリックリンクを表す
- `->`：リンク先を表す

---

## トラブルシューティング

### よくある問題と解決法

#### 1. ダウンロードが途中で止まる

**症状**：ダウンロードが途中で中断される

**原因**：ネットワークエラー、タイムアウト

**解決法**：
```bash
# wgetの再開機能を使う（-cオプション）
wget -c https://...

# または、curlを使う
curl -C - -O https://...

# タイムアウト時間を延長する
wget --timeout=60 -c https://...
```

#### 2. ファイルが見つからない

**症状**：`No such file or directory`

**原因**：パスが間違っている、ファイルが存在しない

**解決法**：
```bash
# 現在のディレクトリを確認
pwd

# ファイルの存在を確認
ls -la /path/to/file

# 絶対パスで指定する
```

#### 3. リネームスクリプトが動かない

**症状**：`python3: command not found` または `Permission denied`

**原因**：Pythonがインストールされていない、またはスクリプトに実行権限がない

**解決法**：
```bash
# Pythonがインストールされているか確認
python3 --version

# 実行権限を付与
chmod +x scripts/rename_chrom.py

# または、python3コマンドで直接実行
python3 scripts/rename_chrom.py ...
```

#### 4. リネームが正しく行われない

**症状**：染色体名がリネームされていない

**原因**：リネーム対応表の形式が間違っている、またはファイル名が一致していない

**解決法**：
```bash
# リネーム対応表を確認
cat chrom_rename.tsv

# FASTAファイルのヘッダーを確認
grep "^>" rios_ancestor.fa | head -5

# 対応表の形式を確認（タブ区切りであることを確認）
```

---

## まとめ

今回は、バックボーンゲノムの取得とデータ準備について解説しました：

1. **Ensemblからのダウンロード**：ゲノムFASTA/GTFファイルの取得方法
2. **染色体リネーム**：架空種のゲノムとして扱いやすくする
3. **シンボリックリンクの作成**：プロジェクトからゲノムファイルを参照

**次回の予告**：
次回は、**トランスクリプト設計とRNA-seqシミュレーション**について詳しく解説します。14遺伝子の選定理由から、polyesterによるリード生成まで、実際のコードを交えながら説明します。

**学習のポイント**：
- ゲノムデータのダウンロードには時間がかかります。途中で止まった場合は、再開機能を使いましょう
- リネームスクリプトは、FASTAとGTFの両方に対応しています
- シンボリックリンクを使うことで、ディスク容量を節約できます

お楽しみに！

---

**参考リンク集**：
- [Ensembl公式サイト](https://www.ensembl.org/)
- [NCBI公式サイト](https://www.ncbi.nlm.nih.gov/)
- [Qiita - ゲノムデータベース関連記事](https://qiita.com/tags/genome)
- [Qiita - Python関連記事](https://qiita.com/tags/python)
- [Qiita - bash関連記事](https://qiita.com/tags/bash)

---
