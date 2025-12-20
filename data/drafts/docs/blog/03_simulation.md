# 第4回：トランスクリプト設計とRNA-seqシミュレーション「polyesterで仮想リードを生成」

## はじめに

前回は、バックボーンゲノムの取得と染色体リネームについて解説しました。今回は、**トランスクリプト設計とRNA-seqシミュレーション**について詳しく説明します。

RNA-seq解析を学ぶには、実際のデータが必要ですが、実験データを取得するのは時間とコストがかかります。そこで、**polyester**というRパッケージを使って、仮想的なRNA-seqリードを生成します。

**この記事で学べること**：
- Rとは何か、なぜバイオインフォマティクスで使われるのか
- polyesterによるRNA-seqリードのシミュレーション方法
- 実験デザインと発現マトリクスの設計方法
- エラーが出たときの対処法

**前提知識**：
- 前回の記事（環境構築）を完了していること
- Rの基本的な知識があると理解しやすいですが、必須ではありません

---

## Rとは？（初心者向け）

### Rとは何か？

**R**は、統計解析やデータ可視化に特化したプログラミング言語です。バイオインフォマティクス分野では、RNA-seq解析の標準ツールとして広く使われています。

**なぜRを使うのか？**

1. **Bioconductor**：バイオインフォマティクス専用のパッケージ集
   - DESeq2、edgeR、polyesterなど、RNA-seq解析に必要なツールが豊富
2. **統計解析に強い**：統計検定、可視化が容易
3. **コミュニティが大きい**：多くのチュートリアルや質問回答がある

**Rの基本操作**：

```r
# Rを起動（ターミナルで）
R

# または、Rスクリプトを実行
Rscript script.R
```

**参考記事**：
- [R入門（初心者向け）](https://qiita.com/tags/r) - QiitaのR関連記事
- [R言語公式サイト](https://www.r-project.org/) - Rの公式サイト
- [Bioconductor公式サイト](https://bioconductor.org/) - Bioconductorの公式サイト

### Rスクリプトの実行方法

**方法1：対話的に実行**：

```bash
# Rを起動
R

# Rのコマンドを入力
> library(polyester)
> # 処理を実行
> q()  # 終了
```

**方法2：スクリプトファイルを実行**：

```bash
# Rスクリプトを実行
Rscript script.R

# または
R --vanilla < script.R
```

**方法3：コマンドライン引数を渡す**：

```bash
# オプションを指定して実行
Rscript script.R --input data.csv --output result.csv
```

**参考記事**：
- [Rスクリプトの実行方法](https://qiita.com/tags/r) - QiitaのRスクリプト記事

---

## 14遺伝子の選定理由

### なぜ14遺伝子だけなのか？

このプロジェクトでは、**14遺伝子のみ**を対象としています。これは、以下の理由からです：

1. **学習目的**：解析手法を理解するため、少数の遺伝子で十分
2. **計算リソース**：メモリとディスク容量を節約
3. **仮説検証**：性差・組織差の主要遺伝子に焦点を当てる

**実際のRNA-seq解析では**：
- 通常、数万〜数十万の遺伝子を対象とする
- しかし、学習目的では少数の遺伝子で十分
- 手法を理解した後、全遺伝子に適用できる

### 選定した遺伝子とその役割

| 遺伝子群 | 遺伝子名 | 役割 | 性差の想定 |
|---------|---------|------|-----------|
| **翼形成** | RiosWingBMP1, RiosWingBMP2 | 翼基部軟骨形成 | オスで高発現 |
| **尾棘形成** | RiosTailSpine1, RiosTailSpine2 | 尾棘形成 | メスで高発現 |
| **毒腺** | RiosToxinA1, RiosToxinA2 | 毒腺分泌タンパク質 | メス特異的 |
| **体軸パターン** | HoxRiosA9, HoxRiosC6 | 体軸パターン形成 | オス：前方シフト、メス：後方シフト |
| **軟骨・骨格** | Sox9Rios, Runx2Rios | 軟骨・骨格形成 | 組織特異的（翼/尾） |
| **色素・模様** | MC1R-Rios, ASIP-Rios, EDNRB-Rios | 色素・模様形成 | 亜種差 |
| **ハウスキーピング** | ACTB_Rios | アクチン（正規化用） | 発現一定 |

### 遺伝子選定の根拠

これらの遺伝子は、性差・形態形成の鍵遺伝子です：

- **Hox遺伝子**：脊椎動物の体軸パターンを決定する遺伝子群
- **BMP/Wnt系**：軟骨・骨格形成に関与
- **Sox9/Runx2**：軟骨・骨格形成のマスター遺伝子
- **MC1R/ASIP/EDNRB**：色素・模様形成に関与（亜種差の検証用）

**参考記事**：
- [遺伝子発現解析の基礎](https://qiita.com/tags/rna-seq) - QiitaのRNA-seq記事
- [Hox遺伝子とは？](https://ja.wikipedia.org/wiki/Hox%E9%81%BA%E4%BC%9D%E5%AD%90) - Wikipedia

<!-- TODO: 画像URLを設定してください: 14遺伝子の選定理由図 -->
<!-- ![14遺伝子の選定理由](画像URL) -->

---

## トランスクリプトFASTAの作成

### 合成配列を使う理由

現実の遺伝子配列を使わず、**ランダム配列（GC含量≈45%）**を生成したのは：

- **「遺伝子機能は配列ではなく発現量で表現される」**という前提
- RNA-seq解析では、配列の生物学的意味より「どの遺伝子がどれだけ発現しているか」が重要
- 配列が現実と違っても、発現マトリクス（TPM/read count）が正しければ解析は成立する

**FASTA形式とは？**

FASTA形式は、DNA/RNA/タンパク質配列を記録する標準的な形式です：

```
>遺伝子名
ATGCGATCGATCGATCG...
ATGCGATCGATCGATCG...
```

- `>`で始まる行：ヘッダー（遺伝子名など）
- その後の行：配列（通常80文字ごとに改行）

**参考記事**：
- [FASTA形式の説明](https://qiita.com/tags/fasta) - QiitaのFASTA記事
- [NCBI FASTA形式の説明](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) - NCBI公式

### トランスクリプトFASTA作成スクリプト

```python
import random

# 14遺伝子の定義
genes = [
    "RiosWingBMP1", "RiosWingBMP2",
    "RiosTailSpine1", "RiosTailSpine2",
    "RiosToxinA1", "RiosToxinA2",
    "HoxRiosA9", "HoxRiosC6",
    "Sox9Rios", "Runx2Rios",
    "MC1R-Rios", "ASIP-Rios", "EDNRB-Rios",
    "ACTB_Rios"
]

# 各遺伝子の長さ（塩基数）
lens = {
    "RiosWingBMP1": 1500, "RiosWingBMP2": 1500,
    "RiosTailSpine1": 1500, "RiosTailSpine2": 1500,
    "RiosToxinA1": 1200, "RiosToxinA2": 1200,
    "HoxRiosA9": 2000, "HoxRiosC6": 2000,
    "Sox9Rios": 1400, "Runx2Rios": 1400,
    "MC1R-Rios": 1100, "ASIP-Rios": 1100, "EDNRB-Rios": 1100,
    "ACTB_Rios": 1800
}

def rand_seq(n, gc=0.45):
    """GC含量を指定してランダム配列を生成"""
    nts = []
    for _ in range(n):
        nts.append(random.choices(
            population=["A", "C", "G", "T"],
            weights=[(1-gc)/2, gc/2, gc/2, (1-gc)/2],
        )[0])
    return "".join(nts)

# FASTAファイルの作成
with open("transcripts.fa", "w") as f:
    for g in genes:
        seq = rand_seq(lens[g])
        f.write(f">{g}\n")
        # 80文字ごとに改行（FASTA形式）
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
```

**スクリプトの実行方法**：

```bash
# Pythonスクリプトを実行
python3 make_transcripts.py

# 生成されたファイルを確認
ls -lh transcripts.fa
head -10 transcripts.fa
```

**ポイント**：
- GC含量を45%に設定（鳥類ゲノムの平均に近い）
- 各遺伝子の長さは、現実の遺伝子長を参考に設定
- FASTA形式で出力（80文字ごとに改行）

**参考記事**：
- [PythonでFASTAファイルを扱う](https://qiita.com/tags/python) - QiitaのPython記事

---

## 実験デザインの設計

### 実験デザインとは？

実験デザインは、**どのサンプルをどの条件で解析するか**を定義する表です。RNA-seq解析では、このデザインに基づいて統計解析を行います。

**実験デザインの重要性**：
- 統計解析の結果に大きく影響する
- 適切なデザインがないと、有意な差を検出できない
- 事前にデザインを決めることが重要

**参考記事**：
- [RNA-seq実験デザインの考え方](https://qiita.com/tags/rna-seq) - Qiitaの実験デザイン記事
- [統計的実験計画法の基礎](https://qiita.com/tags/statistics) - Qiitaの統計記事

### 実験デザインの考え方

19サンプルを設計した理由：

1. **性差検出**：M vs F（例：M_Wing vs F_Wing）
2. **組織差検出**：wing_base vs tail_spine vs toxin_gland
3. **交互作用**：Sex × Tissue（例：オスの翼 vs メスの翼）
4. **亜種差**：volcanicus（翼強化）、sylvestris（毒腺強化）、marina（尾棘強化）

これにより、**「設定通りの性差・組織差が統計的に検出できるか」**を検証できます。

**サンプル数の考え方**：
- 通常、各条件につき3-5サンプル（replicate）が必要
- サンプル数が少ないと、統計的検出力が低下する
- このプロジェクトでは学習目的のため、一部の条件でサンプル数を減らしている

<!-- TODO: 画像URLを設定してください: 実験デザインの可視化図 -->
<!-- ![実験デザインの可視化](画像URL) -->

---

## 発現マトリクスの設計

### 発現マトリクスとは？

発現マトリクスは、**各遺伝子が各サンプルでどれだけ発現しているか**を表す表です。polyesterは、このマトリクスに基づいてリードを生成します。

**発現マトリクスの構造**：

| gene_id | M_Wing_1 | M_Wing_2 | F_Wing_1 | ... |
|---------|----------|----------|----------|-----|
| RiosWingBMP1 | 520 | 500 | 60 | ... |
| RiosWingBMP2 | 480 | 470 | 55 | ... |
| ... | ... | ... | ... | ... |

- **行**：遺伝子
- **列**：サンプル
- **値**：期待されるリード数（TPM相当）

**参考記事**：
- [発現マトリクスの読み方](https://qiita.com/tags/rna-seq) - QiitaのRNA-seq記事

### 発現パターンの設計

設計した発現パターン：

1. **RiosWingBMP群**：
   - オス翼基部：500-520（高発現）
   - メス翼基部：60-70（低発現）
   - オス火山亜種：900（さらに高発現）

2. **RiosTailSpine群**：
   - メス尾棘：780-820（高発現）
   - オス尾棘：110-120（低発現）
   - メス沿岸亜種：880-900（さらに高発現）

3. **RiosToxin群**：
   - メス毒腺：500-540（高発現）
   - オス毒腺：30-35（低発現）
   - メス森林亜種：880-900（さらに高発現）

4. **ACTB_Rios**：
   - 全サンプル：1000（一定、正規化用）

**発現量の単位**：
- **TPM（Transcripts Per Million）**：正規化された発現量
- **Read count**：リード数（生のカウント）
- このプロジェクトでは、TPM相当の値を設定

<!-- TODO: 画像URLを設定してください: 発現マトリクスの可視化 -->
<!-- ![発現マトリクスの可視化](画像URL) -->

---

## polyesterによるリードシミュレーション

### polyesterとは？

**polyester**は、R/Bioconductorパッケージで、**発現マトリクスからRNA-seqリードを生成**するツールです。

**なぜpolyesterを使うのか？**

1. **R/Bioconductor統合**：DESeq2/edgeRと同じ環境で使える
2. **発現量を直接指定できる**：TPMやfold changeを明示的に設定可能
3. **ペアエンド対応**：現実のRNA-seqデータに近い形式
4. **再現性**：乱数種を設定することで、同じ結果を再現できる

**polyesterのインストール**：

```r
# BiocManagerをインストール（初回のみ）
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# polyesterをインストール
BiocManager::install("polyester")
```

**参考記事**：
- [polyester公式ドキュメント](https://bioconductor.org/packages/polyester/) - Bioconductor公式
- [polyesterの使い方](https://qiita.com/tags/polyester) - Qiitaのpolyester記事

### polyester実行スクリプト

```r
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(Biostrings)
  library(polyester)
})

# オプション解析
opt_list <- list(
  make_option("--transcripts", type="character", help="トランスクリプトFASTA"),
  make_option("--design", type="character", help="デザインCSV"),
  make_option("--expr", type="character", help="発現マトリクスCSV"),
  make_option("--out", type="character", default="data/polyester_reads", help="出力ディレクトリ"),
  make_option("--readlen", type="integer", default=150, help="リード長"),
  make_option("--reads_per_sample", type="numeric", default=5e6, help="サンプルあたり総リード数"),
  make_option("--seed", type="integer", default=1234, help="乱数種"),
  make_option("--paired", action="store_true", default=TRUE, help="ペアエンド")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# データ読み込み
design <- fread(opt$design)
expr <- fread(opt$expr)

# 発現マトリクスの準備
expr_mat <- as.matrix(expr[, -1, with=FALSE])
rownames(expr_mat) <- expr$gene_id

# サンプルごとに reads_per_sample に合うようスケーリング
col_sums <- colSums(expr_mat)
scales <- opt$reads_per_sample / col_sums
expr_scaled <- t(t(expr_mat) * scales)
expr_scaled <- round(expr_scaled)

# FASTAファイルの読み込み
fasta <- readDNAStringSet(opt$transcripts)
tx_names <- sub(" .*", "", names(fasta))
names(fasta) <- tx_names

# 一時FASTAファイルを作成（polyesterの要求に合わせる）
temp_fasta_path <- file.path(tempdir(), "clean_transcripts.fa")
writeXStringSet(fasta, temp_fasta_path)

# 発現マトリクスをFASTAの順序に合わせる
readspertx <- expr_scaled[tx_names, , drop=FALSE]

# サンプル数を明示的に指定
num_reps_val <- ncol(readspertx)

# リードシミュレーション実行
simulate_experiment(
  temp_fasta_path,
  reads_per_transcript = readspertx,
  fold_changes = matrix(1, nrow = nrow(readspertx), ncol = num_reps_val),
  outdir = opt$out,
  paired = opt$paired,
  readlen = opt$readlen,
  seed = opt$seed,
  num_reps = num_reps_val
)
```

**スクリプトの説明**：
- `#!/usr/bin/env Rscript`：Rスクリプトであることを示す
- `suppressPackageStartupMessages`：パッケージ読み込み時のメッセージを抑制
- `optparse`：コマンドライン引数を解析
- `data.table`：CSVファイルの読み込み
- `Biostrings`：FASTAファイルの読み込み
- `polyester`：リードシミュレーション

### 実行例

```bash
# conda環境を有効化
conda activate rios-r43

# polyester実行
Rscript simulation/polyester/simulate_polyester.R \
  --transcripts /home/uecha/data/reference/rios/transcripts.fa \
  --design simulation/design/rios_design.csv \
  --expr simulation/polyester/expr_matrix_template.csv \
  --out /home/uecha/scratch/MonsterGenome_Rios/polyester_reads \
  --readlen 150 \
  --reads_per_sample 5e6
```

**実行時間**：
- 19サンプル×14遺伝子×5M reads：数十分〜数時間（マシンスペックによる）

**メモリ使用量**：
- 約4-8GB（サンプル数とリード数による）

### 出力ファイル

polyesterは、以下の形式でリードを出力します：

```
polyester_reads/
├── sample_01_1.fasta  # サンプル1のR1リード（FASTA形式）
├── sample_01_2.fasta  # サンプル1のR2リード（FASTA形式）
├── sample_02_1.fasta
├── sample_02_2.fasta
└── ...
```

**注意**：polyesterは**FASTA形式**で出力しますが、多くのアライメントツール（HISAT2など）は**FASTQ形式**を要求します。後で変換が必要です。

**参考記事**：
- [FASTAとFASTQの違い](https://qiita.com/tags/fastq) - QiitaのFASTQ記事

---

## デバッグ過程の記録

### よくある問題と解決策

#### 1. R環境のバージョン不一致

**問題**：
```
Error: Bioconductor version '3.18' requires R version '4.3'
```

**原因**：Rのバージョンが4.3ではない

**解決策**：
```bash
# R 4.3専用環境を作成（前回の記事を参照）
conda create -y -n rios-r43 -c conda-forge -c bioconda \
  r-base=4.3 bioconductor-polyester=1.38.0

# 環境を有効化
conda activate rios-r43

# Rのバージョンを確認
R --version
```

**参考記事**：
- [Rのバージョン管理](https://qiita.com/tags/r) - QiitaのR記事

#### 2. パッケージ依存の問題

**問題**：
```
Error in library(polyester) : there is no package called 'polyester'
```

**原因**：polyesterがインストールされていない

**解決策**：
```r
# BiocManagerから直接インストール
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("polyester", version = "3.18")

# インストールを確認
library(polyester)
```

**参考記事**：
- [Bioconductorパッケージのインストール方法](https://bioconductor.org/install/) - Bioconductor公式

#### 3. simulate_experimentの引数エラー

**問題**：
```
Error: reads_per_transcript is the wrong length
```

**原因**：FASTAヘッダと発現マトリクスのIDが一致していない

**解決策**：
```r
# FASTAヘッダを確認
fasta <- readDNAStringSet("transcripts.fa")
names(fasta)

# 発現マトリクスのIDを確認
rownames(expr_mat)

# IDを一致させる
tx_names <- sub(" .*", "", names(fasta))
names(fasta) <- tx_names
expr_scaled <- expr_scaled[tx_names, , drop=FALSE]
```

**参考記事**：
- [polyesterのトラブルシューティング](https://bioconductor.org/packages/polyester/) - Bioconductor公式

#### 4. メモリ消費の問題

**問題**：19サンプル×14遺伝子×5M readsでメモリ消費が大きい

**症状**：
- メモリ不足エラー
- 実行が非常に遅い

**解決策**：
```bash
# 軽量版を作成（4サンプル、1M reads、readlen 100bp）
Rscript simulate_polyester.R \
  --reads_per_sample 1e6 \
  --readlen 100 \
  ...

# または、サンプル数を減らす
# デザインCSVを編集して、サンプル数を減らす
```

**参考記事**：
- [Rのメモリ管理](https://qiita.com/tags/r) - QiitaのR記事

---

## 生成されたリードの確認

### リード数の確認

```bash
# 生成されたファイル数の確認
ls /home/uecha/scratch/MonsterGenome_Rios/polyester_reads/*.fasta | wc -l
# 出力例: 38（19サンプル × 2（ペアエンド））

# 最初のサンプルのリード数を確認
grep -c "^>" /home/uecha/scratch/MonsterGenome_Rios/polyester_reads/sample_01_1.fasta
# 出力例: 5000000（5M reads）
```

**コマンドの説明**：
- `ls *.fasta | wc -l`：FASTAファイルの数をカウント
- `grep -c "^>"`：FASTAヘッダー（`>`で始まる行）の数をカウント

### リードの品質確認

```bash
# リード長の確認
head -2 /home/uecha/scratch/MonsterGenome_Rios/polyester_reads/sample_01_1.fasta | tail -1 | wc -c
# 出力例: 151（150bp + 改行）

# リードの内容確認
head -4 /home/uecha/scratch/MonsterGenome_Rios/polyester_reads/sample_01_1.fasta
# 出力例:
# >read_1/1
# ACCCTCAAGTGCAATCTTCTCAAATTGGAGCAAGGTTAAAGGACACATTCTTATGAATGGTTCTGTATATATACGGCACC
# >read_2/1
# ...
```

**確認ポイント**：
- リード数が期待通りか
- リード長が正しいか（150bp）
- ファイルが正しく生成されているか

---

## まとめ

今回は、トランスクリプト設計とRNA-seqシミュレーションについて解説しました：

1. **Rの基本**：Rとは何か、なぜ使うのか
2. **14遺伝子の選定**：性差・形態形成の鍵遺伝子を選定
3. **トランスクリプトFASTA作成**：ランダム配列で合成
4. **実験デザインの設計**：19サンプルで性差・組織差・亜種差を検証
5. **発現マトリクスの設計**：期待される発現パターンを定義
6. **polyesterによるリード生成**：発現マトリクスからリードを生成

**次回の予告**：
次回は、**アライメントとカウント**について詳しく解説します。HISAT2によるアライメントから、トランスクリプトカウントまで、実際のコードを交えながら説明します。

**学習のポイント**：
- Rスクリプトの実行に慣れることが重要です
- エラーが出たときは、エラーメッセージをよく読んで、参考記事を検索してみてください
- メモリ不足の場合は、サンプル数やリード数を減らして試してみてください

お楽しみに！

---

**参考リンク集**：
- [R言語公式サイト](https://www.r-project.org/)
- [Bioconductor公式サイト](https://bioconductor.org/)
- [polyester公式ドキュメント](https://bioconductor.org/packages/polyester/)
- [Qiita - R関連記事](https://qiita.com/tags/r)
- [Qiita - RNA-seq関連記事](https://qiita.com/tags/rna-seq)
- [Qiita - polyester関連記事](https://qiita.com/tags/polyester)

---
