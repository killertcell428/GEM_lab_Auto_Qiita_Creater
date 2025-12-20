# 第7回：可視化と解釈「火山図、ヒートマップ、PCAで見る発現パターン」

## はじめに

前回は、DESeq2による発現解析について解説しました。今回は、**可視化手法と結果の解釈**について詳しく説明します。

RNA-seq解析の結果を理解するには、数値だけでなく、**グラフで可視化**することが重要です。今回は、火山図、ヒートマップ、PCAの3つの可視化手法を紹介します。

---

## パイプライン全体図

この記事は、以下のパイプラインの**ステップ7：可視化**に対応しています：

```
[1] イントロダクション
[2] 環境構築
[3] データ準備
[4] シミュレーション
[5] アライメント・カウント
[6] 発現解析
[7] 可視化 ← 現在ここ
[8] CNV解析
[9] まとめ
```

<!-- TODO: 画像URLを設定してください: パイプライン全体図（ステップ7をハイライト） -->
<!-- ![パイプライン全体図](画像URL) -->

---

## 前回の振り返り

**第6回で学んだこと**：
- 統計検定の基礎（p値、FDR、多重検定補正）
- DESeq2のデザイン式の選択とモデル設計
- 性差・組織差の検定と結果の解釈
- log2FoldChange、pvalue、padjの意味

**今回から始めること**：
- 火山図による発現変動の可視化
- ヒートマップによる発現パターンの可視化
- PCAによるサンプル間の関係性の可視化

---

**この記事で学べること**：
- 可視化の重要性
- ggplot2とは何か、基本的な使い方
- 火山図、ヒートマップ、PCAの読み方
- 結果の解釈方法

**前提知識**：
- 前回の記事（DESeq2解析）を完了していること
- Rの基本的な操作ができること

---

## 可視化の重要性

### なぜ可視化が必要なのか？

**数値だけでは分からないこと**：
- データの分布
- サンプル間の関係性
- 異常値（外れ値）の存在
- パターンや傾向

**可視化の利点**：
- 一目でデータの特徴を把握できる
- 異常値を発見しやすい
- 結果を他人に説明しやすい
- 論文やプレゼンテーションで使える

**参考記事**：
- [データ可視化の重要性](https://qiita.com/tags/visualization) - Qiitaの可視化記事
- [統計グラフの選び方](https://qiita.com/tags/statistics) - Qiitaの統計記事

### ggplot2とは？

**ggplot2**は、Rでグラフを作成するパッケージです。

**ggplot2の特徴**：
- 美しいグラフを簡単に作成できる
- 「グラマー・オブ・グラフィックス」という考え方に基づく
- レイヤー構造でグラフを組み立てる

**ggplot2の基本構造**：
```r
ggplot(data, aes(x=変数1, y=変数2)) +
  geom_point() +  # 点を描画
  labs(title="タイトル") +  # タイトルを設定
  theme_minimal()  # テーマを設定
```

**参考記事**：
- [ggplot2公式ドキュメント](https://ggplot2.tidyverse.org/) - ggplot2公式サイト
- [ggplot2入門](https://qiita.com/tags/ggplot2) - Qiitaのggplot2記事
- [Rでグラフを作成する方法](https://qiita.com/tags/r) - QiitaのR記事

---

## 火山図（Volcano Plot）

### 火山図とは？

火山図は、**発現変動の大きさ（log2FC）と統計的有意性（-log10(pvalue)）**を同時に可視化するグラフです。

**火山図の特徴**：
- 形が火山に似ていることから「火山図」と呼ばれる
- 多くの遺伝子が中央（非有意）に集まり、少数の遺伝子が上部（有意）に分布する

**軸の意味**：
- **X軸**：log2FoldChange（発現変動の大きさ）
- **Y軸**：-log10(pvalue)（統計的有意性）
- **点**：各遺伝子

### 火山図の読み方

**有意領域の判定**：
- **水平線**：pvalue = 0.05（-log10(0.05) ≈ 1.3）
- **垂直線**：log2FC = 0（発現差なし）
- **右上領域**：オスで高発現かつ有意
- **左上領域**：メスで高発現かつ有意

**火山図の見方**：
1. **右上の点**：オスで高発現かつ有意な遺伝子
2. **左上の点**：メスで高発現かつ有意な遺伝子
3. **中央の点**：発現差がない、または非有意な遺伝子

**参考記事**：
- [火山図の読み方](https://qiita.com/tags/volcano-plot) - Qiitaの火山図記事

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

**コードの説明**：
- `ggplot`：グラフの基本設定
- `aes`：軸や色の設定
- `geom_point`：点を描画
- `geom_hline`：水平線を描画（有意水準）
- `geom_vline`：垂直線を描画（発現差なし）
- `ggsave`：グラフをファイルに保存

**参考記事**：
- [ggplot2の基本的な使い方](https://ggplot2.tidyverse.org/) - ggplot2公式

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

**ヒートマップの特徴**：
- 表形式のデータを色で表現
- パターンが一目で分かる

**軸の意味**：
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

**コードの説明**：
- `melt`：データをlong形式に変換（ggplot2で扱いやすい形式）
- `scale`：z-score正規化（平均0、標準偏差1に変換）
- `geom_tile`：タイル（四角形）を描画
- `scale_fill_gradient2`：3色のグラデーション（青→白→赤）

**z-score正規化とは？**
- 各遺伝子の発現量を、平均0、標準偏差1に変換
- 遺伝子間で発現量のスケールが異なる場合に有効

**参考記事**：
- [z-score正規化とは？](https://qiita.com/tags/normalization) - Qiitaの正規化記事

### ヒートマップの読み方

**カラースケール**：
- **青**：その遺伝子の平均発現量より低い
- **白**：その遺伝子の平均発現量
- **赤**：その遺伝子の平均発現量より高い

**パターンの解釈**：
- **RiosWingBMP1/2**：オス翼サンプル（M_Wing_*）で赤色（高発現）
- **RiosToxinA1/2**：メス毒腺サンプル（F_Toxin_*）で赤色（高発現）
- **RiosTailSpine1/2**：メス尾棘サンプル（F_Tail_*）で赤色（高発現）

**参考記事**：
- [ヒートマップの読み方](https://qiita.com/tags/heatmap) - Qiitaのヒートマップ記事

<!-- TODO: 画像URLを設定してください: ヒートマップ -->
<!-- ![ヒートマップ: 全遺伝子の発現パターン](画像URL) -->

---

## PCA（Principal Component Analysis）

### PCAとは？

**PCA（主成分分析）**は、**サンプル間の関係性**を可視化する手法です。高次元のデータを2-3次元に圧縮し、サンプルがどのようにグループ化されるかを確認できます。

**なぜPCAが必要なのか？**
- 14遺伝子×19サンプルのデータは14次元
- 14次元のデータを直接可視化することはできない
- PCAにより、2-3次元に圧縮して可視化できる

**PCAの原理**：
- データの変動が大きい方向（主成分）を見つける
- 第1主成分（PC1）：最も大きな変動を説明する軸
- 第2主成分（PC2）：2番目に大きな変動を説明する軸

**参考記事**：
- [PCAとは何か？](https://qiita.com/tags/pca) - QiitaのPCA記事
- [主成分分析の基礎](https://qiita.com/tags/statistics) - Qiitaの統計記事

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

**コードの説明**：
- `prcomp`：PCAを計算する関数
- `t(expr_matrix)`：行列を転置（サンプル×遺伝子に変換）
- `scale. = TRUE`：標準化（平均0、標準偏差1に変換）
- `pca_data$x[, 1:2]`：PC1とPC2の値を取得

### PCAの読み方

**主成分（PC）の意味**：
- **PC1**：最も大きな変動を説明する軸
- **PC2**：2番目に大きな変動を説明する軸
- **寄与率**：各PCが全体の変動の何%を説明するか

**解釈**：
- **近いサンプル**：発現パターンが似ている
- **遠いサンプル**：発現パターンが異なる
- **グループ化**：組織や性別でグループ化されることが多い

**参考記事**：
- [PCAの読み方](https://qiita.com/tags/pca) - QiitaのPCA記事

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

## トラブルシューティング

### よくある問題と解決法

#### 1. ggplot2がインストールされていない

**症状**：`Error in library(ggplot2) : there is no package called 'ggplot2'`

**解決策**：
```r
# ggplot2をインストール
install.packages("ggplot2")

# または、conda環境でインストール
conda install -n rios-env r-ggplot2
```

#### 2. グラフが保存されない

**症状**：`ggsave`でエラーが出る

**解決策**：
```r
# ディレクトリが存在するか確認
dir.create("analysis/out/plots", recursive=TRUE)

# 絶対パスで保存
ggsave("/full/path/to/plot.png", plot, ...)
```

#### 3. グラフが表示されない

**症状**：グラフが表示されない、または空白

**解決策**：
```r
# グラフを直接表示
print(volcano_plot)

# または、対話的に実行
R
> source("visualize.R")
> print(volcano_plot)
```

---

## まとめ

今回は、可視化手法と結果の解釈について解説しました：

1. **可視化の重要性**：数値だけでは分からないパターンを可視化
2. **ggplot2の基本**：Rでグラフを作成する方法
3. **火山図**：発現変動の大きさと統計的有意性を可視化
4. **ヒートマップ**：遺伝子×サンプルの発現パターンを可視化
5. **PCA**：サンプル間の関係性を可視化
6. **主要所見**：設計通りの性差・組織差が確認された

**次回の予告**：
次回は、**CNV解析とスライド用データ作成**について詳しく解説します。CNV領域の定義から、期待値と検出値の比較まで、実際のコードを交えながら説明します。

**学習のポイント**：
- 可視化は、データを理解する上で重要です
- ggplot2は、美しいグラフを簡単に作成できます
- グラフの読み方を理解することで、結果をより深く理解できます

お楽しみに！

---

**参考リンク集**：
- [ggplot2公式ドキュメント](https://ggplot2.tidyverse.org/)
- [Qiita - ggplot2関連記事](https://qiita.com/tags/ggplot2)
- [Qiita - 可視化関連記事](https://qiita.com/tags/visualization)
- [Qiita - PCA関連記事](https://qiita.com/tags/pca)

---
