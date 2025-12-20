# 第5回：アライメントとカウント「HISAT2とトランスクリプトカウント」

## はじめに

前回は、polyesterによるRNA-seqリードのシミュレーションについて解説しました。今回は、生成したリードを**アライメント（マッピング）**し、**カウント（リード数を数える）**する過程について詳しく説明します。

RNA-seq解析では、リードをリファレンスゲノム（またはトランスクリプト）にマッピングし、各遺伝子にマップされたリード数を数えることで、発現量を定量化します。

---

## パイプライン全体図

この記事は、以下のパイプラインの**ステップ5：アライメント・カウント**に対応しています：

```
[1] イントロダクション
[2] 環境構築
[3] データ準備
[4] シミュレーション
[5] アライメント・カウント ← 現在ここ
[6] 発現解析
[7] 可視化
[8] CNV解析
[9] まとめ
```

<!-- TODO: 画像URLを設定してください: パイプライン全体図（ステップ5をハイライト） -->
<!-- ![パイプライン全体図](画像URL) -->

---

## 前回の振り返り

**第4回で学んだこと**：
- Rの基本操作とBioconductorの使い方
- 14遺伝子の選定理由とトランスクリプトFASTA作成
- 実験デザイン（19サンプル）の設計と発現マトリクスの作成
- polyesterによるRNA-seqリード生成（約380M reads）

**今回から始めること**：
- FASTA→FASTQ変換
- HISAT2トランスクリプトインデックスの作成
- アライメントの実行とBAMファイルの生成
- カスタムカウント関数によるリードカウント

---

**この記事で学べること**：
- アライメントとは何か、なぜ必要なのか
- HISAT2とは何か、どのように使うのか
- BAMファイルとは何か
- リードカウントの方法

**前提知識**：
- 前回の記事（シミュレーション）を完了していること
- ターミナルの基本的な操作ができること

---

## アライメントとは？（初心者向け）

### アライメント（マッピング）とは？

**アライメント（マッピング）**は、**リード配列をリファレンス配列（ゲノムやトランスクリプト）に「どこに位置するか」を特定する**処理です。

**なぜアライメントが必要なのか？**
- RNA-seqリードは、どの遺伝子から来たのかが分からない
- アライメントすることで、各リードがどの遺伝子に対応するかを特定できる
- 遺伝子ごとのリード数を数えることで、発現量を定量化できる

**アライメントのイメージ**：
```
リード配列: ATGCGATCGATCG...
リファレンス: ...ATGCGATCGATCGATCGATCG...
              ↑ ここにマップされる
```

**参考記事**：
- [アライメント（マッピング）とは？](https://qiita.com/tags/alignment) - Qiitaのアライメント記事
- [RNA-seq解析の基礎](https://qiita.com/tags/rna-seq) - QiitaのRNA-seq記事

### HISAT2とは？

**HISAT2**は、RNA-seqリードをゲノムやトランスクリプトにマッピングするツールです。

**HISAT2の特徴**：
- 高速：大量のリードを短時間でマッピングできる
- メモリ効率が良い：インデックスが軽量
- ペアエンド対応：2本のリード（R1/R2）を同時に処理できる

**HISAT2のインストール**：
```bash
# conda環境を有効化
conda activate rios-env

# HISAT2がインストールされているか確認
which hisat2
# 出力例: /home/uecha/miniconda3/envs/rios-env/bin/hisat2

# バージョンを確認
hisat2 --version
```

**参考記事**：
- [HISAT2公式ドキュメント](https://daehwankimlab.github.io/hisat2/) - HISAT2公式サイト
- [HISAT2の使い方](https://qiita.com/tags/hisat2) - QiitaのHISAT2記事

---

## アライメント戦略の選択

### ゲノムアライメント vs トランスクリプトアライメント

RNA-seq解析では、2つのアライメント戦略があります：

1. **ゲノムアライメント**
   - リードをゲノム配列（`rios_genome.fa`）にマッピング
   - スプライスジャンクションを考慮する必要がある
   - ツール：STAR, HISAT2（ゲノムインデックス）

2. **トランスクリプトアライメント**
   - リードをトランスクリプト配列（`transcripts.fa`）にマッピング
   - スプライスジャンクションを考慮する必要がない
   - ツール：HISAT2（トランスクリプトインデックス）、Salmon、kallisto

**スプライスジャンクションとは？**
- 真核生物の遺伝子は、エクソンとイントロンで構成される
- スプライシングにより、イントロンが除去され、エクソンが結合される
- ゲノムアライメントでは、このスプライシングを考慮する必要がある

**参考記事**：
- [スプライシングとは？](https://ja.wikipedia.org/wiki/%E3%82%B9%E3%83%97%E3%83%A9%E3%82%A4%E3%82%B9) - Wikipedia
- [RNA-seqアライメント戦略の比較](https://qiita.com/tags/rna-seq) - QiitaのRNA-seq記事

### このプロジェクトでの選択

このプロジェクトでは、**トランスクリプトアライメント**を選択しました。理由は：

- polyesterで生成したリードは、`transcripts.fa`から直接生成されている
- スプライスジャンクションがないため、解析がシンプル
- トランスクリプトインデックスは軽量で、メモリ消費が少ない

**参考記事**：
- [トランスクリプトアライメント vs ゲノムアライメント](https://qiita.com/tags/rna-seq) - QiitaのRNA-seq記事

<!-- TODO: 画像URLを設定してください: アライメント戦略の比較図 -->
<!-- ![アライメント戦略の比較図](画像URL) -->

---

## FASTA→FASTQ変換

### なぜ変換が必要なのか？

polyesterは**FASTA形式**でリードを出力しますが、多くのアライメントツール（HISAT2など）は**FASTQ形式**を要求します。

**FASTA形式とは？**
- DNA/RNA配列を記録する形式
- ヘッダー行（`>`で始まる）と配列行で構成

**FASTQ形式とは？**
- DNA/RNA配列と品質スコアを記録する形式
- 4行で1リードを表現（ID行、配列行、`+`行、品質スコア行）

**FASTA形式の例**：
```
>read_1/1
ACCCTCAAGTGCAATCTTCTCAAATTGGAGCAAGGTTAAAGGACACATTCTTATGAATGGTTCTGTATATATACGGCACC
>read_2/1
...
```

**FASTQ形式の例**：
```
@read_1/1
ACCCTCAAGTGCAATCTTCTCAAATTGGAGCAAGGTTAAAGGACACATTCTTATGAATGGTTCTGTATATATACGGCACC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read_2/1
...
```

**品質スコアとは？**
- 各塩基の読み取り精度を表す
- Phred score：エラー確率の対数変換
- このプロジェクトでは、シミュレーションデータなので品質スコアは固定値

**参考記事**：
- [FASTAとFASTQの違い](https://qiita.com/tags/fastq) - QiitaのFASTQ記事
- [品質スコアとは？](https://qiita.com/tags/quality-score) - Qiitaの品質スコア記事

### 変換スクリプトの実装

```python
#!/usr/bin/env python3
"""FASTA形式のリードをFASTQ形式に変換（品質スコアを固定値で埋める）"""

import sys
from pathlib import Path

def fasta_to_fastq(fasta_path, fastq_path, quality='I'):
    """FASTAをFASTQに変換（品質スコアは全て同じ値）"""
    with open(fasta_path) as fin, open(fastq_path, 'w') as fout:
        seq_id = None
        seq = []
        for line in fin:
            line = line.rstrip()
            if line.startswith('>'):
                # 前の配列を出力
                if seq_id and seq:
                    seq_str = ''.join(seq)
                    fout.write(f"@{seq_id}\n{seq_str}\n+\n{quality * len(seq_str)}\n")
                # 新しい配列開始
                seq_id = line[1:].split()[0]  # >以降、最初の空白まで
                seq = []
            else:
                seq.append(line)
        # 最後の配列
        if seq_id and seq:
            seq_str = ''.join(seq)
            fout.write(f"@{seq_id}\n{seq_str}\n+\n{quality * len(seq_str)}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 convert_fasta_to_fastq.py input.fasta output.fastq")
        sys.exit(1)
    
    fasta_path = Path(sys.argv[1])
    fastq_path = Path(sys.argv[2])
    
    fasta_to_fastq(fasta_path, fastq_path)
    print(f"Converted: {fasta_path} -> {fastq_path}")
```

**スクリプトの実行方法**：
```bash
# Pythonスクリプトを実行
python3 convert_fasta_to_fastq.py \
  sample_01_1.fasta \
  sample_01_1.fastq

# 変換結果を確認
head -4 sample_01_1.fastq
```

**ポイント**：
- 品質スコアは固定値（'I' = Phred score 40）を使用
- シミュレーションデータなので、品質スコアの意味は重要ではない

**参考記事**：
- [Pythonでファイルを変換する方法](https://qiita.com/tags/python) - QiitaのPython記事

---

## トランスクリプトインデックスの作成

### インデックスとは？

アライメントツールは、リファレンス配列を**事前にインデックス化**する必要があります。インデックスを作成することで、リードを高速にマッピングできます。

**インデックスのイメージ**：
- 本の索引のようなもの
- 「この配列はどこにあるか」を素早く検索できる

**インデックス作成の時間**：
- トランスクリプト（14遺伝子）：数秒〜数分
- ゲノム全体：数時間〜数日

### HISAT2トランスクリプトインデックスの作成

```bash
#!/bin/bash
# HISAT2トランスクリプトインデックスの作成

TRANSCRIPTS_FA="/home/uecha/data/reference/rios/transcripts.fa"
INDEX_DIR="/home/uecha/data/databases/rios/hisat2_transcript_index"
INDEX_PREFIX="${INDEX_DIR}/transcripts"

mkdir -p "${INDEX_DIR}"

echo "[info] Building HISAT2 transcript index..."
hisat2-build -p 4 "${TRANSCRIPTS_FA}" "${INDEX_PREFIX}"

echo "[done] Index created at ${INDEX_PREFIX}"
```

**コマンドの説明**：
- `hisat2-build`：HISAT2のインデックス作成コマンド
- `-p 4`：4スレッドを使用（並列処理）
- 最初の引数：リファレンスFASTAファイル
- 2番目の引数：インデックスファイルのプレフィックス（拡張子なし）

**実行例**：
```bash
# スクリプトを保存（scripts/build_transcript_index.sh）
# 実行権限を付与
chmod +x scripts/build_transcript_index.sh

# スクリプトを実行
bash scripts/build_transcript_index.sh
```

**生成されるファイル**：
```
hisat2_transcript_index/
├── transcripts.1.ht2
├── transcripts.2.ht2
├── transcripts.3.ht2
├── transcripts.4.ht2
├── transcripts.5.ht2
├── transcripts.6.ht2
├── transcripts.7.ht2
└── transcripts.8.ht2
```

**参考記事**：
- [HISAT2インデックスの作成方法](https://daehwankimlab.github.io/hisat2/) - HISAT2公式

---

## HISAT2アライメントの実行

### アライメントコマンド

```bash
# ペアエンドリードのアライメント
hisat2 -p 4 \
  -x /home/uecha/data/databases/rios/hisat2_transcript_index/transcripts \
  -1 sample_01_1.fastq \
  -2 sample_01_2.fastq \
  | samtools sort -@ 4 -o M_Wing_1.bam -
```

**コマンドの説明**：
- `hisat2`：HISAT2のアライメントコマンド
- `-p 4`：4スレッドを使用
- `-x`：インデックスファイルのパス（拡張子なし）
- `-1`, `-2`：ペアエンドリードのR1/R2ファイル
- `|`：パイプライン（前のコマンドの出力を次のコマンドの入力に渡す）
- `samtools sort`：BAMファイルをソート（後続の解析に必要）
- `-o M_Wing_1.bam`：出力ファイル名
- `-`：標準入力から読み込む（パイプライン経由）

**パイプライン（`|`）とは？**
- 前のコマンドの出力を、次のコマンドの入力に渡す
- 例：`hisat2 ... | samtools sort ...` → HISAT2の出力をsamtoolsに渡す

**参考記事**：
- [パイプライン（`|`）の使い方](https://qiita.com/tags/pipe) - Qiitaのパイプライン記事

### BAMファイルとは？

**BAMファイル**は、アライメント結果を記録する形式です。

**BAMファイルの特徴**：
- バイナリ形式（人間が読めない）
- SAMファイルを圧縮した形式
- ファイルサイズが小さい

**BAMファイルの確認方法**：
```bash
# samtoolsでBAMファイルを確認
samtools view M_Wing_1.bam | head -5

# BAMファイルの統計情報を確認
samtools flagstat M_Wing_1.bam
```

**samtoolsとは？**
- BAM/SAMファイルを操作するツール
- ソート、インデックス作成、統計情報の取得などができる

**参考記事**：
- [BAM/SAMファイルとは？](https://qiita.com/tags/bam) - QiitaのBAM記事
- [samtools公式ドキュメント](https://www.htslib.org/) - samtools公式サイト

### アライメント率の確認

```bash
# アライメント統計の確認
samtools flagstat M_Wing_1.bam

# 出力例:
# 10000000 + 0 in total (QC-passed reads + QC-failed reads)
# 9999999 + 0 primary mapped (100.00% : N/A)
# 0 + 0 secondary
# 0 + 0 supplementary
# 0 + 0 duplicates
# 9999999 + 0 mapped (100.00% : N/A)
```

**出力の見方**：
- `in total`：総リード数
- `primary mapped`：マッピングされたリード数
- `mapped`：マッピング率（%）

**期待される結果**：
- マッピング率：99.99%以上（トランスクリプトから生成したリードなので、ほぼ100%マップされる）

**参考記事**：
- [samtools flagstatの読み方](https://qiita.com/tags/samtools) - Qiitaのsamtools記事

---

## リードカウント

### リードカウントとは？

**リードカウント**は、**各遺伝子にマップされたリード数を数える**処理です。

**なぜリードカウントが必要なのか？**
- リード数 = 遺伝子の発現量の指標
- リード数が多い = その遺伝子が多く発現している

**カウントマトリクスとは？**
- 遺伝子×サンプルのリード数を表す表
- 行：遺伝子
- 列：サンプル
- 値：リード数

### featureCountsの限界

通常、RNA-seq解析では`featureCounts`を使ってリードをカウントします。しかし、このプロジェクトでは**トランスクリプトアライメント**を使用しているため、`featureCounts`は使えません。

**理由**：
- `featureCounts`は、ゲノム座標に基づいてリードをカウントする
- トランスクリプトアライメントのBAMには、ゲノム座標がない
- 代わりに、トランスクリプトIDがリファレンス名として記録されている

**参考記事**：
- [featureCountsとは？](https://qiita.com/tags/featurecounts) - QiitaのfeatureCounts記事

### カスタムカウント関数の実装

```r
# トランスクリプトベースのカウント関数
count_transcripts <- function(bam_files) {
  library(Rsamtools)
  library(data.table)
  
  # 各BAMファイルからトランスクリプトIDを抽出
  counts_list <- list()
  
  for (bam_file in bam_files) {
    sample_name <- gsub("\\.bam$", "", basename(bam_file))
    
    # BAMファイルを読み込み
    bam <- scanBam(bam_file, param=ScanBamParam(
      what=c("rname", "qname")
    ))[[1]]
    
    # リファレンス名（トランスクリプトID）をカウント
    tx_counts <- table(bam$rname)
    
    counts_list[[sample_name]] <- as.numeric(tx_counts)
    names(counts_list[[sample_name]]) <- names(tx_counts)
  }
  
  # データフレームに変換
  all_tx <- unique(unlist(lapply(counts_list, names)))
  count_matrix <- matrix(0, nrow=length(all_tx), ncol=length(counts_list))
  rownames(count_matrix) <- all_tx
  colnames(count_matrix) <- names(counts_list)
  
  for (sample_name in names(counts_list)) {
    count_matrix[names(counts_list[[sample_name]]), sample_name] <- 
      counts_list[[sample_name]]
  }
  
  # 出力形式に整形
  dt_out <- data.table(
    Geneid = rownames(count_matrix),
    Chr = "transcript",
    Start = 1,
    End = 1,
    Strand = "+",
    Length = 1,
    count_matrix
  )
  
  return(dt_out)
}
```

**スクリプトの説明**：
- `Rsamtools`：RでBAMファイルを読み込むパッケージ
- `scanBam`：BAMファイルを読み込む関数
- `table`：各トランスクリプトIDの出現回数をカウント
- `data.table`：データフレームを効率的に操作するパッケージ

**実行方法**：
```r
# Rスクリプトを実行
Rscript count_transcripts.R

# または、Rで対話的に実行
R
> source("count_transcripts.R")
> bam_files <- list.files("analysis/out", pattern="\\.bam$", full.names=TRUE)
> counts <- count_transcripts(bam_files)
> fwrite(counts, "analysis/out/gene_counts.txt", sep="\t")
```

**ポイント**：
- `Rsamtools`パッケージでBAMファイルを読み込み
- `rname`（リファレンス名）がトランスクリプトID
- 各トランスクリプトIDの出現回数をカウント

**参考記事**：
- [Rsamtools公式ドキュメント](https://bioconductor.org/packages/Rsamtools/) - Bioconductor公式
- [RでBAMファイルを扱う](https://qiita.com/tags/r) - QiitaのR記事

<!-- TODO: 画像URLを設定してください: パイプライン実行フローチャート -->
<!-- ![パイプライン実行フローチャート](画像URL) -->

---

## 実行結果の確認

### カウントマトリクスの確認

```bash
# カウントファイルの確認
head -5 analysis/out/gene_counts.txt

# 出力例:
# Geneid	Chr	Start	End	Strand	Length	M_Wing_1	M_Wing_2	...
# ACTB_Rios	transcript	1	1	+	1	1416997	1405119	...
# ASIP-Rios	transcript	1	1	+	1	113706	119396	...
# EDNRB-Rios	transcript	1	1	+	1	141861	146760	...
```

**カウントマトリクスの見方**：
- `Geneid`：遺伝子名
- `M_Wing_1`：サンプル名
- 数値：リード数

### アライメント率の統計

```r
# 全サンプルのアライメント率を確認
library(data.table)

bam_files <- list.files("analysis/out", pattern="\\.bam$", full.names=TRUE)
alignment_rates <- sapply(bam_files, function(bam) {
  stats <- system(paste("samtools flagstat", bam), intern=TRUE)
  # マッピング率を抽出
  mapped_line <- grep("mapped", stats, value=TRUE)
  # パースして返す
  # ...
})

mean(alignment_rates)  # 出力例: 99.99
```

---

## よくある問題と解決策

### 1. アライメント率が0%

**症状**：すべてのリードがマッピングされない

**原因**：
- ゲノムインデックスに対してトランスクリプトリードをアライメントしていた
- インデックスが正しく作成されていない
- リードファイルのパスが間違っている

**解決策**：
```bash
# インデックスファイルの存在を確認
ls -lh /home/uecha/data/databases/rios/hisat2_transcript_index/

# インデックスを再作成
hisat2-build -p 4 transcripts.fa index_prefix

# リードファイルのパスを確認
ls -lh sample_01_1.fastq
```

**参考記事**：
- [HISAT2のトラブルシューティング](https://daehwankimlab.github.io/hisat2/) - HISAT2公式

### 2. featureCountsでカウント0%

**症状**：`featureCounts`でカウントが0になる

**原因**：
- トランスクリプトアライメントのBAMにはゲノム座標がない
- `featureCounts`はゲノム座標を要求する

**解決策**：
- カスタムカウント関数を使用
- BAMファイルから直接トランスクリプトIDを抽出

**参考記事**：
- [トランスクリプトアライメントのカウント方法](https://qiita.com/tags/rna-seq) - QiitaのRNA-seq記事

### 3. メモリ不足

**症状**：メモリ不足エラー、実行が非常に遅い

**原因**：
- 大量のBAMファイルを同時に読み込んでいる
- インデックス作成時にメモリが不足

**解決策**：
```bash
# サンプルごとに順次処理
for bam in *.bam; do
  # 処理を実行
done

# インデックス作成時はメモリ制限を設定
ulimit -v 8000000  # 8GBに制限
```

**参考記事**：
- [メモリ管理の方法](https://qiita.com/tags/memory) - Qiitaのメモリ記事

### 4. パイプラインエラー

**症状**：`|`を使ったコマンドが動かない

**原因**：PowerShellでは`|`の動作が異なる場合がある

**解決策**：
```bash
# 中間ファイルを使う
hisat2 ... > temp.sam
samtools sort temp.sam -o output.bam

# または、bashで実行
bash -c "hisat2 ... | samtools sort ..."
```

---

## まとめ

今回は、アライメントとカウントについて解説しました：

1. **アライメントとは何か**：リードをリファレンスにマッピングする処理
2. **HISAT2の使い方**：インデックス作成とアライメント実行
3. **BAMファイルとは何か**：アライメント結果の記録形式
4. **リードカウント**：カスタム関数でトランスクリプトIDをカウント
5. **トラブルシューティング**：よくある問題と解決策

**次回の予告**：
次回は、**発現解析と統計検定**について詳しく解説します。DESeq2の使い方から、モデル設計の考え方まで、実際のコードを交えながら説明します。

**学習のポイント**：
- アライメントは、RNA-seq解析の中核となる処理です
- インデックスを作成することで、高速にアライメントできます
- エラーが出たときは、インデックスファイルやリードファイルのパスを確認してください

お楽しみに！

---

**参考リンク集**：
- [HISAT2公式ドキュメント](https://daehwankimlab.github.io/hisat2/)
- [samtools公式ドキュメント](https://www.htslib.org/)
- [Rsamtools公式ドキュメント](https://bioconductor.org/packages/Rsamtools/)
- [Qiita - HISAT2関連記事](https://qiita.com/tags/hisat2)
- [Qiita - samtools関連記事](https://qiita.com/tags/samtools)
- [Qiita - アライメント関連記事](https://qiita.com/tags/alignment)

---
