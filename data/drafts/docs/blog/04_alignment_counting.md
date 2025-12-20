# 第5回：アライメントとカウント「HISAT2とトランスクリプトカウント」

## はじめに

前回は、polyesterによるRNA-seqリードのシミュレーションについて解説しました。今回は、生成したリードを**アライメント（マッピング）**し、**カウント（リード数を数える）**する過程について詳しく説明します。

RNA-seq解析では、リードをリファレンスゲノム（またはトランスクリプト）にマッピングし、各遺伝子にマップされたリード数を数えることで、発現量を定量化します。

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

### このプロジェクトでの選択

このプロジェクトでは、**トランスクリプトアライメント**を選択しました。理由は：

- polyesterで生成したリードは、`transcripts.fa`から直接生成されている
- スプライスジャンクションがないため、解析がシンプル
- トランスクリプトインデックスは軽量で、メモリ消費が少ない

<!-- TODO: 画像URLを設定してください: アライメント戦略の比較図 -->
<!-- ![アライメント戦略の比較図](画像URL) -->

---

## FASTA→FASTQ変換

### なぜ変換が必要なのか？

polyesterは**FASTA形式**でリードを出力しますが、多くのアライメントツール（HISAT2など）は**FASTQ形式**を要求します。

**FASTA形式**：
```
>read_1/1
ACCCTCAAGTGCAATCTTCTCAAATTGGAGCAAGGTTAAAGGACACATTCTTATGAATGGTTCTGTATATATACGGCACC
>read_2/1
...
```

**FASTQ形式**：
```
@read_1/1
ACCCTCAAGTGCAATCTTCTCAAATTGGAGCAAGGTTAAAGGACACATTCTTATGAATGGTTCTGTATATATACGGCACC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read_2/1
...
```

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

**ポイント**：
- 品質スコアは固定値（'I' = Phred score 40）を使用
- シミュレーションデータなので、品質スコアの意味は重要ではない

---

## トランスクリプトインデックスの作成

### インデックスとは？

アライメントツールは、リファレンス配列を**事前にインデックス化**する必要があります。インデックスを作成することで、リードを高速にマッピングできます。

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

**実行例**：
```bash
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

**パラメータの説明**：
- `-p 4`：4スレッドを使用
- `-x`：インデックスファイルのパス（拡張子なし）
- `-1`, `-2`：ペアエンドリードのR1/R2ファイル
- `samtools sort`：BAMファイルをソート（後続の解析に必要）

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

**期待される結果**：
- マッピング率：99.99%以上（トランスクリプトから生成したリードなので、ほぼ100%マップされる）

---

## リードカウント

### featureCountsの限界

通常、RNA-seq解析では`featureCounts`を使ってリードをカウントします。しかし、このプロジェクトでは**トランスクリプトアライメント**を使用しているため、`featureCounts`は使えません。

**理由**：
- `featureCounts`は、ゲノム座標に基づいてリードをカウントする
- トランスクリプトアライメントのBAMには、ゲノム座標がない
- 代わりに、トランスクリプトIDがリファレンス名として記録されている

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

**ポイント**：
- `Rsamtools`パッケージでBAMファイルを読み込み
- `rname`（リファレンス名）がトランスクリプトID
- 各トランスクリプトIDの出現回数をカウント

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

**原因**：
- ゲノムインデックスに対してトランスクリプトリードをアライメントしていた
- インデックスが正しく作成されていない

**解決策**：
- トランスクリプトインデックスを作成
- インデックスファイルの存在を確認

### 2. featureCountsでカウント0%

**原因**：
- トランスクリプトアライメントのBAMにはゲノム座標がない
- `featureCounts`はゲノム座標を要求する

**解決策**：
- カスタムカウント関数を使用
- BAMファイルから直接トランスクリプトIDを抽出

### 3. メモリ不足

**原因**：
- 大量のBAMファイルを同時に読み込んでいる
- インデックス作成時にメモリが不足

**解決策**：
- サンプルごとに順次処理
- インデックス作成時はメモリ制限を設定

---

## まとめ

今回は、アライメントとカウントについて解説しました：

1. **アライメント戦略の選択**：トランスクリプトアライメントを選択
2. **FASTA→FASTQ変換**：polyesterの出力をアライメントツール用に変換
3. **トランスクリプトインデックスの作成**：HISAT2でインデックスを作成
4. **アライメントの実行**：HISAT2でリードをマッピング
5. **リードカウント**：カスタム関数でトランスクリプトIDをカウント

次回は、**発現解析と統計検定**について詳しく解説します。DESeq2の使い方から、モデル設計の考え方まで、実際のコードを交えながら説明します。

お楽しみに！

---
