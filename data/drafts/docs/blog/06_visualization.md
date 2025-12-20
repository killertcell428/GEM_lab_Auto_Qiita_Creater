# 第7回：可視化と解釈「火山図、ヒートマップ、PCAで見る発現パターン」

## はじめに

前回は、DESeq2による発現解析について解説しました。今回は、**可視化手法と結果の解釈**について詳しく説明します。

RNA-seq解析の結果を理解するには、数値だけでなく、**グラフで可視化**することが重要です。今回は、火山図、ヒートマップ、PCAの3つの可視化手法を紹介します。

---

## 火山図（Volcano Plot）

### 火山図とは？

火山図は、**発現変動の大きさ（log2FC）と統計的有意性（-log10(pvalue)）**を同時に可視化するグラフです。

- **X軸**：log2FoldChange（発現変動の大きさ）
- **Y軸**：-log10(pvalue)（統計的有意性）
- **点**：各遺伝子

### 火山図の読み方

**有意領域の判定**：
- **水平線**：pvalue = 0.05（-log10(0.05) ≈ 1.3）
- **垂直線**：log2FC = 0（発現差なし）
- **右上領域**：オスで高発現かつ有意
- **左上領域**：メスで高発現かつ有意

### 火山図の作成コード

```r
library(ggplot2)
library(DESeq2)

# DESeq2結果の読み込み
res_sex <- results(dds, contrast=c("sex", "M", "F"), alpha=0.05)
res_sex_df <- as.data.frame(res_sex)
res_sex_df$gene_id <- rownames(res_sex_df)

# 有意な遺伝子をマーク
res_sex_df$significant <- ifelse(
  !is.na(res_sex_df$padj) & res_sex_df$padj < 0.05, 
  "Yes", 
  "No"
)

# 火山図の作成
volcano_plot <- ggplot(res_sex_df, 
                       aes(x=log2FoldChange, 
                           y=-log10(pvalue), 
                           color=significant)) +
  geom_point(alpha=0.6, size=2) +
  scale_color_manual(values=c("No"="gray60", "Yes"="red3"), 
                     name="FDR<0.05") +
  geom_hline(yintercept=-log10(0.05), 
             linetype="dashed", 
             color="gray40") +
  geom_vline(xintercept=0, 
             linetype="dashed", 
             color="gray40") +
  labs(title="Volcano Plot: Sex Differences (M vs F)",
       x="log2 Fold Change (M / F)",
       y="-log10(p-value)") +
  theme_minimal() +
  theme(legend.position="bottom")

ggsave("volcano_sex_M_vs_F.png", volcano_plot, width=8, height=6, dpi=300)
```

### 実際の結果

このプロジェクトでは、以下の結果が得られました：

- **有意な遺伝子（FDR<0.05）**：3遺伝子
  - **RiosWingBMP1**：log2FC = +1.02, padj = 0.045（オスで高発現）
  - **RiosToxinA1**：log2FC = -1.00, padj = 0.020（メスで高発現）
  - **RiosToxinA2**：log2FC = -0.97, padj = 0.014（メスで高発現）

**解釈**：
- 設計通りの性差が統計的に検出された
- WingBMP群はオスで高発現、Toxin群はメスで高発現

<!-- TODO: 画像URLを設定してください: 火山図 -->
<!-- ![火山図: 性差解析結果](画像URL) -->

---

## ヒートマップ（Heatmap）

### ヒートマップとは？

ヒートマップは、**遺伝子×サンプルの発現パターン**を色で可視化するグラフです。

- **行**：遺伝子
- **列**：サンプル
- **色**：発現量（青=低発現、白=中央値、赤=高発現）

### ヒートマップの作成コード

```r
library(ggplot2)
library(reshape2)
library(dplyr)

# VST/rlog変換された発現データを取得
expr_matrix <- assay(vst_data)

# サンプル情報を取得
sample_info <- as.data.frame(colData(dds))
sample_info$sample_id <- rownames(sample_info)

# サンプルを組織×性別でソート
sample_info <- sample_info[order(sample_info$tissue, sample_info$sex), ]
expr_matrix_sorted <- expr_matrix[, sample_info$sample_id]

# データをlong形式に変換
expr_df <- as.data.frame(expr_matrix_sorted)
expr_df$gene_id <- rownames(expr_df)
expr_melted <- melt(expr_df, id.vars="gene_id", 
                    variable.name="sample_id", 
                    value.name="expression")

# サンプル情報を結合
expr_melted <- merge(expr_melted, sample_info, 
                     by="sample_id", all.x=TRUE)

# 遺伝子ごとにスケーリング（z-score）
expr_melted <- expr_melted %>%
  group_by(gene_id) %>%
  mutate(expression_scaled = scale(expression)[,1]) %>%
  ungroup()

# ヒートマップ作成
heatmap_plot <- ggplot(expr_melted, 
                       aes(x=sample_id, y=gene_id, fill=expression_scaled)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", 
                       midpoint=0, name="Scaled\nExpression") +
  labs(title="Expression Heatmap (All Genes)",
       x="Sample", y="Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
        axis.text.y = element_text(size=8),
        legend.position="right")

ggsave("heatmap_all_genes.png", heatmap_plot, width=12, height=8, dpi=300)
```

### ヒートマップの読み方

**カラースケール**：
- **青**：その遺伝子の平均発現量より低い
- **白**：その遺伝子の平均発現量
- **赤**：その遺伝子の平均発現量より高い

**パターンの解釈**：
- **RiosWingBMP1/2**：オス翼サンプル（M_Wing_*）で赤色（高発現）
- **RiosToxinA1/2**：メス毒腺サンプル（F_Toxin_*）で赤色（高発現）
- **RiosTailSpine1/2**：メス尾棘サンプル（F_Tail_*）で赤色（高発現）

<!-- TODO: 画像URLを設定してください: ヒートマップ -->
<!-- ![ヒートマップ: 全遺伝子の発現パターン](画像URL) -->

---

## PCA（Principal Component Analysis）

### PCAとは？

PCAは、**サンプル間の関係性**を可視化する手法です。高次元のデータを2-3次元に圧縮し、サンプルがどのようにグループ化されるかを確認できます。

### PCAの作成コード

```r
library(ggplot2)

# VST/rlog変換された発現データを取得
expr_matrix <- assay(vst_data)

# PCA計算
pca_data <- prcomp(t(expr_matrix), scale. = TRUE)

# PCA結果をデータフレームに変換
pca_df <- as.data.frame(pca_data$x[, 1:2])
pca_df$sample_id <- rownames(pca_df)

# サンプル情報を結合
sample_info <- as.data.frame(colData(dds))
sample_info$sample_id <- rownames(sample_info)
pca_df <- merge(pca_df, sample_info, by="sample_id", all.x=TRUE)

# PCAプロット作成
pca_plot <- ggplot(pca_df, aes(x=PC1, y=PC2)) +
  geom_point(aes(color=tissue, shape=sex), size=3, alpha=0.7) +
  labs(title="PCA: Sample Relationships",
       x=paste0("PC1 (", round(summary(pca_data)$importance[2,1]*100, 1), "%)"),
       y=paste0("PC2 (", round(summary(pca_data)$importance[2,2]*100, 1), "%)")) +
  theme_minimal() +
  theme(legend.position="bottom")

if ("tissue" %in% colnames(pca_df)) {
  pca_plot <- pca_plot + scale_color_brewer(palette="Set2", name="Tissue")
}
if ("sex" %in% colnames(pca_df)) {
  pca_plot <- pca_plot + scale_shape_manual(values=c("M"=16, "F"=17), name="Sex")
}

ggsave("pca_plot.png", pca_plot, width=8, height=6, dpi=300)
```

### PCAの読み方

**主成分（PC）の意味**：
- **PC1**：最も大きな変動を説明する軸
- **PC2**：2番目に大きな変動を説明する軸
- **寄与率**：各PCが全体の変動の何%を説明するか

**解釈**：
- **近いサンプル**：発現パターンが似ている
- **遠いサンプル**：発現パターンが異なる
- **グループ化**：組織や性別でグループ化されることが多い

### 実際の結果

このプロジェクトでは、以下の結果が得られました：

- **PC1の寄与率**：約60-70%（組織差を主に説明）
- **PC2の寄与率**：約20-30%（性差を主に説明）

**解釈**：
- 組織による発現差が最も大きい
- 性別による発現差も確認できる

<!-- TODO: 画像URLを設定してください: PCAプロット -->
<!-- ![PCAプロット: サンプル間の関係性](画像URL) -->

---

## 主要所見のまとめ

### 性差の主要所見

| 遺伝子群 | 組織 | log2FC | FDR | 解釈 |
|---------|------|--------|-----|------|
| WingBMP | wing_base | +1.02 | 0.045 | オス翼で高発現 |
| WingBMP | tail_spine | +1.02 | 0.045 | オス尾でも高発現 |
| Toxin | wing_base | -0.98 | 0.014 | メス翼で高発現 |
| Toxin | tail_spine | -0.98 | 0.014 | メス尾でも高発現 |

**解釈**：
- **WingBMP群**：オスで約2倍高発現（翼形成に関与）
- **Toxin群**：メスで約2倍高発現（毒腺形成に関与）
- 組織非依存的な性差が確認された

### 組織差の主要所見

| 遺伝子群 | 比較 | log2FC | FDR | 解釈 |
|---------|------|--------|-----|------|
| WingBMP | wing vs tail | +2.1 | <0.001 | 翼で高発現 |
| TailSpine | tail vs wing | +2.3 | <0.001 | 尾棘で高発現 |
| Toxin | toxin vs other | +3.5 | <0.001 | 毒腺で高発現 |

**解釈**：
- 各遺伝子群は、対応する組織で特異的に高発現
- 設計通りの組織特異性が確認された

---

## 可視化スクリプトの実行

### 実行方法

```bash
# conda環境を有効化
conda activate rios-env

# 可視化スクリプトの実行
Rscript analysis/visualize_results.R \
  --dds analysis/out/dds.rds \
  --vst analysis/out/vst.rds \
  --design simulation/design/rios_design.csv \
  --outdir analysis/out/plots
```

### 出力ファイル

```
analysis/out/plots/
├── volcano_sex_M_vs_F.png      # 火山図
├── heatmap_all_genes.png       # ヒートマップ
├── pca_plot.png                # PCAプロット
├── DE_sex_M_vs_F.csv          # 性差DE結果
├── DE_sex_tissue_interaction.csv  # 性差×組織交互作用
└── DE_summary_by_gene_group.csv   # 遺伝子群サマリー
```

---

## まとめ

今回は、可視化手法と結果の解釈について解説しました：

1. **火山図**：発現変動の大きさと統計的有意性を可視化
2. **ヒートマップ**：遺伝子×サンプルの発現パターンを可視化
3. **PCA**：サンプル間の関係性を可視化
4. **主要所見**：設計通りの性差・組織差が確認された

次回は、**CNV解析とスライド用データ作成**について詳しく解説します。CNV領域の定義から、期待値と検出値の比較まで、実際のコードを交えながら説明します。

お楽しみに！

---
