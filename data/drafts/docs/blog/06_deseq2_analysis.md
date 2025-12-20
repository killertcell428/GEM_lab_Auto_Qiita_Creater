# 第6回：発現解析と統計検定「DESeq2で性差と組織差を検出」

## はじめに

前回は、アライメントとカウントについて解説しました。今回は、**DESeq2による発現解析**について詳しく説明します。

DESeq2は、RNA-seqデータの統計解析を行うRパッケージです。サンプル間の発現差を検定し、どの遺伝子が有意に変動しているかを判定します。

---

## パイプライン全体図

この記事は、以下のパイプラインの**ステップ6：発現解析**に対応しています：

```
[1] イントロダクション
[2] 環境構築
[3] データ準備
[4] シミュレーション
[5] アライメント・カウント
[6] 発現解析 ← 現在ここ
[7] 可視化
[8] CNV解析
[9] まとめ
```

<!-- TODO: 画像URLを設定してください: パイプライン全体図（ステップ6をハイライト） -->
<!-- ![パイプライン全体図](画像URL) -->

---

## 前回の振り返り

**第5回で学んだこと**：
- アライメント（マッピング）とは何か、なぜ必要なのか
- HISAT2トランスクリプトインデックスの作成方法
- アライメントの実行とBAMファイルの生成
- カスタムカウント関数によるリードカウント（カウントマトリクスの作成）

**今回から始めること**：
- DESeq2による統計解析の実行
- デザイン式の選択とモデル設計
- 性差・組織差の検定と結果の解釈

---

**この記事で学べること**：
- DESeq2とは何か、なぜ使うのか
- 統計検定の基礎（p値、FDRとは何か）
- デザイン式の考え方
- 結果の解釈方法

**前提知識**：
- 前回の記事（アライメントとカウント）を完了していること
- Rの基本的な操作ができること
- 統計の基礎知識があると理解しやすいですが、必須ではありません

---

## 統計検定の基礎（初心者向け）

### 統計検定とは？

**統計検定**は、**「2つのグループに差があるか」を統計的に判定する**手法です。

**なぜ統計検定が必要なのか？**
- データには偶然の変動（ノイズ）が含まれる
- 「見た目の差」が本当に意味のある差なのかを判定する必要がある
- 統計検定により、偶然の差と意味のある差を区別できる

**統計検定の流れ**：
1. **仮説を立てる**：例「オスとメスで発現量に差がある」
2. **検定を実行**：データから統計量を計算
3. **p値を計算**：偶然でこのような差が出る確率
4. **判定**：p値が小さい（通常<0.05）なら「有意な差がある」

**参考記事**：
- [統計検定の基礎](https://qiita.com/tags/statistics) - Qiitaの統計記事
- [p値とは何か？](https://qiita.com/tags/p-value) - Qiitaのp値記事

### p値とは？

**p値**は、**「偶然でこのような結果が出る確率」**を表します。

**p値の解釈**：
- **p値 < 0.05**：偶然でこのような差が出る確率が5%未満 → 「有意な差がある」
- **p値 ≥ 0.05**：偶然でこのような差が出る確率が5%以上 → 「有意な差がない」

**なぜ0.05なのか？**
- 統計学の慣習で、5%を「有意水準」として使うことが多い
- 厳しくする場合は0.01、緩くする場合は0.10を使うこともある

**p値の注意点**：
- p値が小さい = 「差がある」ことを示すが、「差の大きさ」は示さない
- p値が大きい = 「差がない」ことを示すが、「本当に差がない」とは言えない

**参考記事**：
- [p値の正しい理解](https://qiita.com/tags/p-value) - Qiitaのp値記事
- [統計的有意性とは？](https://qiita.com/tags/statistics) - Qiitaの統計記事

### 多重検定補正（FDR）とは？

**多重検定補正**は、**複数の遺伝子を同時に検定する際、偽陽性を補正する**手法です。

**なぜ多重検定補正が必要なのか？**
- 14遺伝子を検定する場合、p値<0.05で判定すると、偶然で1つは有意になる可能性がある
- 1000遺伝子を検定する場合、p値<0.05で判定すると、偶然で50個は有意になる可能性がある
- 多重検定補正により、偽陽性を減らせる

**FDR（False Discovery Rate）とは？**
- 多重検定補正の手法の1つ
- 「有意と判定した遺伝子のうち、どれくらいが偽陽性か」を制御する
- **padj（adjusted p-value）**：FDR補正後のp値

**FDRの解釈**：
- **padj < 0.05**：FDRが5%未満 → 「有意な差がある」
- **padj ≥ 0.05**：FDRが5%以上 → 「有意な差がない」

**参考記事**：
- [多重検定補正とは？](https://qiita.com/tags/multiple-testing) - Qiitaの多重検定記事
- [FDRとは何か？](https://qiita.com/tags/fdr) - QiitaのFDR記事

---

## DESeq2とは？

### DESeq2の役割

**DESeq2**は、RNA-seqデータの統計解析を行うRパッケージです。

DESeq2は、以下の処理を行います：

1. **正規化**：サンプル間のリード総数の差を補正（size factor）
2. **分散推定**：遺伝子ごとの発現変動をモデル化
3. **検定**：コントラスト（例：M_Wing vs F_Wing）で t検定相当
4. **多重検定補正**：FDR（False Discovery Rate）で調整

**なぜDESeq2を使うのか？**

RNA-seqデータの特徴：

- **カウントデータ**：離散値（0, 1, 2, ...）
- **分散が大きい**：発現量が高い遺伝子ほど分散が大きい
- **サンプル数が少ない**：通常3-5サンプル程度

DESeq2は、これらの特徴を考慮した統計モデルを使用します。

**参考記事**：
- [DESeq2公式ドキュメント](https://bioconductor.org/packages/DESeq2/) - Bioconductor公式
- [DESeq2の使い方](https://qiita.com/tags/deseq2) - QiitaのDESeq2記事
- [RNA-seq解析の統計手法](https://qiita.com/tags/rna-seq) - QiitaのRNA-seq記事

<!-- TODO: 画像URLを設定してください: DESeq2の処理フロー図 -->
<!-- ![DESeq2の処理フロー図](画像URL) -->

---

## デザイン式の選択

### デザイン式とは？

デザイン式は、**どの因子（性別、組織、亜種など）を考慮するか**を定義する式です。

**デザイン式の例**：
- `~ sex`：性別のみを考慮
- `~ tissue`：組織のみを考慮
- `~ sex + tissue`：性別と組織の両方を考慮
- `~ sex + tissue + sex:tissue`：性別、組織、交互作用を考慮

**交互作用とは？**
- 性別による発現差が組織によって異なること
- 例：オスの翼では高発現だが、オスの尾では低発現

**参考記事**：
- [実験デザインの考え方](https://qiita.com/tags/experimental-design) - Qiitaの実験デザイン記事
- [交互作用とは？](https://qiita.com/tags/interaction) - Qiitaの交互作用記事

### このプロジェクトでのデザイン式

このプロジェクトでは、以下のデザイン式を使用しました：

```r
design_formula <- ~ sex + tissue + sex:tissue
```

**意味**：
- `sex`：性別の主効果（M vs F）
- `tissue`：組織の主効果（wing_base vs tail_spine vs ...）
- `sex:tissue`：性別×組織の交互作用（例：オスの翼 vs メスの翼）

**デザイン式の選択方法**：
1. データの構造を確認（どの因子があるか）
2. 仮説を確認（どの因子の効果を検証したいか）
3. サンプル数を確認（複雑なデザインには多くのサンプルが必要）

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

**参考記事**：
- [フルランク問題の対処法](https://qiita.com/tags/deseq2) - QiitaのDESeq2記事

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

**コードの説明**：
- `DESeqDataSetFromMatrix`：カウントマトリクスとメタデータからDESeq2オブジェクトを作成
- `factor`：性別や組織を因子（カテゴリ変数）として扱う
- `rowSums(counts(dds)) >= 10`：全サンプルで合計10リード以上の遺伝子のみを残す
- `DESeq`：DESeq2のメイン解析関数

**実行時間**：
- 14遺伝子×19サンプル：数秒〜数分

**参考記事**：
- [DESeq2の基本的な使い方](https://bioconductor.org/packages/DESeq2/) - Bioconductor公式

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

**出力の見方**：
- `LFC > 0 (up)`：オスで高発現した遺伝子数
- `LFC < 0 (down)`：メスで高発現した遺伝子数
- `outliers`：外れ値として検出された遺伝子数

---

## 性差解析の結果

### 性差（M vs F）のコントラスト

**コントラストとは？**
- 2つのグループを比較する設定
- 例：`contrast=c("sex", "M", "F")` → オス（M）とメス（F）を比較

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

**交互作用の検証方法**：
- 各組織ごとに性差を計算
- 組織間で性差の大きさが異なるか確認

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

**変換前の問題**：
- 発現量が高い遺伝子ほど分散が大きい
- 可視化やクラスタリングが困難

**変換後の利点**：
- 分散が安定化される
- 可視化やクラスタリングが容易になる

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

**参考記事**：
- [VSTとrlogの違い](https://qiita.com/tags/deseq2) - QiitaのDESeq2記事

---

## 統計検定結果の解釈

### log2FoldChange（log2FC）

**log2FoldChangeとは？**
- 発現量の比を対数（底2）で表した値
- オスとメスの発現量の比を表す

**なぜ対数変換するのか？**
- 発現量の比は、対数変換することで正規分布に近づく
- 統計検定が容易になる

**意味**：
- 正の値：比較群（M）で高発現
- 負の値：比較群（M）で低発現
- 絶対値が大きいほど、発現差が大きい

**例**：
- log2FC = +1.0 → 2倍高発現（2^1.0 = 2.0）
- log2FC = -1.0 → 2倍低発現（2^-1.0 = 0.5）
- log2FC = +2.0 → 4倍高発現（2^2.0 = 4.0）

**参考記事**：
- [Fold Changeとは？](https://qiita.com/tags/fold-change) - QiitaのFold Change記事
- [対数変換とは？](https://qiita.com/tags/log-transformation) - Qiitaの対数変換記事

### pvalueとpadj（FDR）

**pvalue**：
- その遺伝子が偶然でこのような発現差を示す確率
- 通常、pvalue < 0.05を有意とみなす

**padj（FDR）**：
- 多重検定補正後のp値
- 複数の遺伝子を同時に検定するため、偽陽性を補正
- 通常、padj < 0.05を有意とみなす

**pvalueとpadjの違い**：
- **pvalue**：その遺伝子単独で検定した場合のp値
- **padj**：複数の遺伝子を同時に検定した場合の補正後p値
- **padjの方が厳しい基準**：padj < 0.05の方が、pvalue < 0.05より厳しい

**参考記事**：
- [p値とFDRの違い](https://qiita.com/tags/statistics) - Qiitaの統計記事

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

**判定基準**：
- **padj < 0.05**：統計的有意性
- **|log2FC| > 1.0**：生物学的に意味のある差（2倍以上）

**参考記事**：
- [有意な遺伝子の選び方](https://qiita.com/tags/deseq2) - QiitaのDESeq2記事

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

## トラブルシューティング

### よくある問題と解決法

#### 1. フルランクエラー

**症状**：`model matrix is not full rank`

**原因**：サンプル数が少ない、またはすべての組み合わせが存在しない

**解決策**：
```r
# より単純なデザインを試す
design_formula <- ~ sex + tissue  # 交互作用を除外
# または
design_formula <- ~ tissue  # 性別も除外
```

#### 2. メモリ不足

**症状**：メモリ不足エラー

**原因**：遺伝子数やサンプル数が多い

**解決策**：
```r
# 発現量が低い遺伝子を除外
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
```

#### 3. 結果が空

**症状**：`results()`で結果が取得できない

**原因**：コントラストの指定が間違っている

**解決策**：
```r
# 因子のレベルを確認
levels(colData(dds)$sex)
levels(colData(dds)$tissue)

# 正しいレベル名でコントラストを指定
res <- results(dds, contrast=c("sex", "M", "F"))
```

---

## まとめ

今回は、DESeq2による発現解析について解説しました：

1. **統計検定の基礎**：p値、FDRとは何か
2. **DESeq2の使い方**：デザイン式の選択と実行
3. **結果の解釈**：log2FC、pvalue、padjの意味
4. **性差・組織差の解析**：コントラストによる検定
5. **VST vs rlog**：遺伝子数に応じた選択

**次回の予告**：
次回は、**可視化と解釈**について詳しく解説します。火山図、ヒートマップ、PCAによる発現パターンの可視化と、結果の解釈方法を説明します。

**学習のポイント**：
- 統計検定の基礎を理解することが重要です
- p値とFDRの違いを理解してください
- デザイン式は、データの構造に合わせて選択してください

お楽しみに！

---

**参考リンク集**：
- [DESeq2公式ドキュメント](https://bioconductor.org/packages/DESeq2/)
- [DESeq2論文](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
- [Qiita - DESeq2関連記事](https://qiita.com/tags/deseq2)
- [Qiita - 統計検定関連記事](https://qiita.com/tags/statistics)
- [Qiita - RNA-seq関連記事](https://qiita.com/tags/rna-seq)

---
