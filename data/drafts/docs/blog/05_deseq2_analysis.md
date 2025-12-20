# 第6回：発現解析と統計検定「DESeq2で性差と組織差を検出」

## はじめに

前回は、アライメントとカウントについて解説しました。今回は、**DESeq2による発現解析**について詳しく説明します。

DESeq2は、RNA-seqデータの統計解析を行うRパッケージです。サンプル間の発現差を検定し、どの遺伝子が有意に変動しているかを判定します。

---

## DESeq2とは？

### DESeq2の役割

DESeq2は、以下の処理を行います：

1. **正規化**：サンプル間のリード総数の差を補正（size factor）
2. **分散推定**：遺伝子ごとの発現変動をモデル化
3. **検定**：コントラスト（例：M_Wing vs F_Wing）で t検定相当
4. **多重検定補正**：FDR（False Discovery Rate）で調整

### なぜDESeq2を使うのか？

RNA-seqデータの特徴：

- **カウントデータ**：離散値（0, 1, 2, ...）
- **分散が大きい**：発現量が高い遺伝子ほど分散が大きい
- **サンプル数が少ない**：通常3-5サンプル程度

DESeq2は、これらの特徴を考慮した統計モデルを使用します。

<!-- TODO: 画像URLを設定してください: DESeq2の処理フロー図 -->
<!-- ![DESeq2の処理フロー図](画像URL) -->

---

## デザイン式の選択

### デザイン式とは？

デザイン式は、**どの因子（性別、組織、亜種など）を考慮するか**を定義する式です。

### このプロジェクトでのデザイン式

このプロジェクトでは、以下のデザイン式を使用しました：

```r
design_formula <- ~ sex + tissue + sex:tissue
```

**意味**：
- `sex`：性別の主効果（M vs F）
- `tissue`：組織の主効果（wing_base vs tail_spine vs ...）
- `sex:tissue`：性別×組織の交互作用（例：オスの翼 vs メスの翼）

### モデル行列のフルランク問題

**フルランクとは？**

モデル行列がフルランクであるとは、すべてのパラメータが独立に推定できることを意味します。

**問題が発生する場合**：

- 亜種サンプルが各1つしかない場合
- すべての組み合わせが存在しない場合（例：オスの毒腺サンプルがない）

**解決策**：

```r
# より単純なデザインを試す
if (grepl("not full rank", error_message)) {
  design_formula <- ~ sex + tissue  # 交互作用を除外
  # または
  design_formula <- ~ tissue  # 性別も除外
}
```

<!-- TODO: 画像URLを設定してください: デザイン式の選択フローチャート -->
<!-- ![デザイン式の選択フローチャート](画像URL) -->

---

## DESeq2の実行

### 基本的な実行コード

```r
library(DESeq2)
library(data.table)

# カウントファイルの読み込み
counts <- fread("analysis/out/gene_counts.txt")
gene_ids <- counts$Geneid
count_matrix <- as.matrix(counts[, -c(1:6)])

# メタデータの準備
col_data <- fread("simulation/design/rios_design.csv")
col_data <- col_data[match(colnames(count_matrix), col_data$sample_id), ]
col_data$sex <- factor(col_data$sex)
col_data$tissue <- factor(col_data$tissue)

# DESeqDataSetの作成
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ sex + tissue + sex:tissue
)

# フィルタリング（発現量が低い遺伝子を除外）
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# DESeq2の実行
dds <- DESeq(dds)
```

### 実行結果の確認

```r
# 結果のサマリー
summary(dds)

# 出力例:
# out of 14 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1, 7.1%
# LFC < 0 (down)     : 2, 14.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
```

---

## 性差解析の結果

### 性差（M vs F）のコントラスト

```r
# 性差の結果を取得
res_sex <- results(dds, contrast=c("sex", "M", "F"), alpha=0.05)

# 結果をデータフレームに変換
res_sex_df <- as.data.frame(res_sex)
res_sex_df$gene_id <- rownames(res_sex_df)

# 有意な遺伝子を抽出
sig_genes <- res_sex_df[!is.na(res_sex_df$padj) & res_sex_df$padj < 0.05, ]
```

### 主要な結果

| 遺伝子 | log2FC | pvalue | padj | 解釈 |
|--------|--------|--------|------|------|
| RiosWingBMP1 | +1.02 | 0.0096 | 0.045 | オスで高発現 |
| RiosToxinA1 | -1.00 | 0.0028 | 0.020 | メスで高発現 |
| RiosToxinA2 | -0.97 | 0.0010 | 0.014 | メスで高発現 |

**解釈**：
- **RiosWingBMP1**：オスで約2倍（2^1.02 ≈ 2.0）高発現
- **RiosToxinA1/A2**：メスで約2倍（2^-1.0 ≈ 0.5）高発現

---

## 組織差解析の結果

### 組織差（wing_base vs tail_spine）のコントラスト

```r
# 組織差の結果を取得
res_tissue <- results(dds, contrast=c("tissue", "wing_base", "tail_spine"), alpha=0.05)

# 結果の確認
head(res_tissue[order(res_tissue$padj), ])
```

### 主要な結果

| 遺伝子 | log2FC | padj | 解釈 |
|--------|--------|------|------|
| RiosWingBMP1 | +2.1 | <0.001 | 翼で高発現 |
| RiosTailSpine1 | -2.3 | <0.001 | 尾棘で高発現 |

---

## 性差×組織交互作用の解析

### 交互作用とは？

交互作用とは、**性別による発現差が組織によって異なる**ことを意味します。

例：
- オスの翼：RiosWingBMP1が高発現
- メスの翼：RiosWingBMP1が低発現
- オスの尾：RiosWingBMP1が低発現
- メスの尾：RiosWingBMP1が低発現

### 各組織での性差を計算

```r
# 各組織での性差を計算
tissues <- unique(colData(dds)$tissue)
sex_tissue_results <- list()

for (tiss in tissues) {
  # その組織のサンプルのみを抽出
  dds_tiss <- dds[, colData(dds)$tissue == tiss]
  
  # 性別が両方存在するかチェック
  if (length(unique(colData(dds_tiss)$sex)) >= 2) {
    res_st <- results(dds_tiss, contrast=c("sex", "M", "F"), alpha=0.05)
    res_st_df <- as.data.frame(res_st)
    res_st_df$tissue <- tiss
    sex_tissue_results[[tiss]] <- res_st_df
  }
}
```

### 主要な結果

| 遺伝子群 | 組織 | log2FC | padj | 解釈 |
|---------|------|--------|------|------|
| WingBMP | wing_base | +1.02 | 0.045 | オス翼で高発現 |
| WingBMP | tail_spine | +1.02 | 0.045 | オス尾でも高発現（組織非依存） |
| Toxin | wing_base | -0.98 | 0.014 | メス翼で高発現 |
| Toxin | tail_spine | -0.98 | 0.014 | メス尾でも高発現（組織非依存） |

---

## VST vs rlogの選択

### なぜ変換が必要なのか？

RNA-seqデータは、発現量が高い遺伝子ほど分散が大きいという特徴があります。可視化やクラスタリングのため、**分散を安定化**する変換が必要です。

### VST（Variance Stabilizing Transformation）

```r
# VST変換
vst_data <- varianceStabilizingTransformation(dds, blind=FALSE)
```

**特徴**：
- 高速
- 遺伝子数が多い場合（≥1000）に適している

### rlog（Regularized log transformation）

```r
# rlog変換
rlog_data <- rlog(dds, blind=FALSE)
```

**特徴**：
- 遺伝子数が少ない場合（<1000）に適している
- より正確だが、計算時間がかかる

### このプロジェクトでの選択

このプロジェクトでは、**遺伝子数が14と少ない**ため、`rlog`を使用しました：

```r
n_genes <- nrow(counts)
if (n_genes < 1000) {
  message("[info] Gene count is low. Using rlog transformation.")
  vst_data <- rlog(dds, blind=FALSE)
} else {
  vst_data <- varianceStabilizingTransformation(dds, blind=FALSE)
}
```

---

## 統計検定結果の解釈

### log2FoldChange（log2FC）

**意味**：
- 正の値：比較群（M）で高発現
- 負の値：比較群（M）で低発現
- 絶対値が大きいほど、発現差が大きい

**例**：
- log2FC = +1.0 → 2倍高発現（2^1.0 = 2.0）
- log2FC = -1.0 → 2倍低発現（2^-1.0 = 0.5）
- log2FC = +2.0 → 4倍高発現（2^2.0 = 4.0）

### pvalueとpadj（FDR）

**pvalue**：
- その遺伝子が偶然でこのような発現差を示す確率
- 通常、pvalue < 0.05を有意とみなす

**padj（FDR）**：
- 多重検定補正後のp値
- 複数の遺伝子を同時に検定するため、偽陽性を補正
- 通常、padj < 0.05を有意とみなす

### 有意な遺伝子の判定

```r
# 有意な遺伝子を抽出
sig_genes <- res_sex_df[
  !is.na(res_sex_df$padj) & 
  res_sex_df$padj < 0.05 & 
  abs(res_sex_df$log2FoldChange) > 1.0,  # log2FC > 1.0（2倍以上の差）
]

# 統計サマリー
n_sig <- nrow(sig_genes)
n_up_m <- sum(sig_genes$log2FoldChange > 0)
n_up_f <- sum(sig_genes$log2FoldChange < 0)

message("Significant genes (FDR<0.05, |log2FC|>1.0): ", n_sig)
message("  - Up in M: ", n_up_m)
message("  - Up in F: ", n_up_f)
```

---

## 主要遺伝子群のサマリー

### 遺伝子群ごとの発現変動

| gene_group | tissue | avg_log2FC | min_FDR | n_genes |
|-----------|--------|------------|---------|---------|
| WingBMP | wing_base | 1.02 | 0.045 | 2 |
| WingBMP | tail_spine | 1.02 | 0.045 | 2 |
| Toxin | wing_base | -0.98 | 0.014 | 2 |
| Toxin | tail_spine | -0.98 | 0.014 | 2 |
| TailSpine | wing_base | -0.61 | 0.064 | 2 |
| TailSpine | tail_spine | -0.61 | 0.064 | 2 |

**解釈**：
- **WingBMP群**：オスで約2倍高発現（log2FC = +1.02）
- **Toxin群**：メスで約2倍高発現（log2FC = -0.98）
- **TailSpine群**：メスで約1.5倍高発現（log2FC = -0.61）

---

## まとめ

今回は、DESeq2による発現解析について解説しました：

1. **デザイン式の選択**：sex + tissue + sex:tissue
2. **モデル行列のフルランク問題**：サンプル数が少ない場合の対処
3. **性差・組織差の解析**：コントラストによる検定
4. **VST vs rlog**：遺伝子数に応じた選択
5. **結果の解釈**：log2FC、pvalue、padjの意味

次回は、**可視化と解釈**について詳しく解説します。火山図、ヒートマップ、PCAによる発現パターンの可視化と、結果の解釈方法を説明します。

お楽しみに！

---
