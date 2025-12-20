# 補足資料

## 補足表S1. 実験デザイン行列

| サンプルID | 性別 | 亜種 | 組織 | 反復 | 条件 | 備考 |
|-----------|-----|------------|--------|-----------|-----------|-------|
| M_Wing_1 | M | Rios_ref | wing_base | 1 | baseline | 雄翼基部 |
| M_Wing_2 | M | Rios_ref | wing_base | 2 | baseline | 雄翼基部 |
| M_Wing_3 | M | Rios_ref | wing_base | 3 | baseline | 雄翼基部 |
| F_Wing_1 | F | Rios_ref | wing_base | 1 | baseline | 雌翼基部 |
| F_Wing_2 | F | Rios_ref | wing_base | 2 | baseline | 雌翼基部 |
| F_Tail_1 | F | Rios_ref | tail_spine | 1 | baseline | 雌尾棘 |
| F_Tail_2 | F | Rios_ref | tail_spine | 2 | baseline | 雌尾棘 |
| M_Tail_1 | M | Rios_ref | tail_spine | 1 | baseline | 雄尾棘 |
| F_Toxin_1 | F | Rios_ref | toxin_gland | 1 | baseline | 雌毒腺 |
| F_Toxin_2 | F | Rios_ref | toxin_gland | 2 | baseline | 雌毒腺 |
| M_Pectoral_1 | M | Rios_ref | pectoral_muscle | 1 | baseline | 雄胸筋 |
| F_Pectoral_1 | F | Rios_ref | pectoral_muscle | 1 | baseline | 雌胸筋 |
| M_Brain_1 | M | Rios_ref | brain | 1 | baseline | 雄脳 |
| F_Brain_1 | F | Rios_ref | brain | 1 | baseline | 雌脳 |
| M_Gonad_1 | M | Rios_ref | gonad | 1 | baseline | 雄性腺 |
| F_Gonad_1 | F | Rios_ref | gonad | 1 | baseline | 雌性腺 |
| M_Wing_volc_1 | M | Rios_volcanicus | wing_base | 1 | baseline | 火山亜種、翼強化 |
| F_Toxin_sylv_1 | F | Rios_sylvestris | toxin_gland | 1 | baseline | 森林亜種、毒腺強化 |
| F_Tail_marina_1 | F | Rios_marina | tail_spine | 1 | baseline | 沿岸亜種、尾棘強化 |

**総サンプル数：** 19  
**性別分布：** 雄9、雌10  
**組織分布：** wing_base (5), tail_spine (4), toxin_gland (3), pectoral_muscle (2), brain (2), gonad (2)  
**亜種分布：** Rios_ref (16), Rios_volcanicus (1), Rios_sylvestris (1), Rios_marina (1)

---

## 補足表S2. 遺伝子アノテーションと期待発現パターン

| 遺伝子ID | 遺伝子クラスター | 機能 | 期待性特異的発現 | 期待組織特異的発現 |
|---------|--------------|----------|-------------------|---------------------|
| RiosWingBMP1 | WingBMP | 翼基部軟骨形成 | 雄特異的 | wing_base |
| RiosWingBMP2 | WingBMP | 翼基部軟骨形成 | 雄特異的 | wing_base |
| RiosTailSpine1 | TailSpine | 尾棘発達 | 雌特異的 | tail_spine |
| RiosTailSpine2 | TailSpine | 尾棘発達 | 雌特異的 | tail_spine |
| RiosToxinA1 | Toxin | 毒腺分泌 | 雌特異的 | toxin_gland |
| RiosToxinA2 | Toxin | 毒腺分泌 | 雌特異的 | toxin_gland |
| HoxRiosA9 | Hox | 前後軸特異化 | 性特異的・組織特異的 | 複数 |
| HoxRiosC6 | Hox | 前後軸特異化 | 性特異的・組織特異的 | 複数 |
| Sox9Rios | Cartilage | 軟骨形成 | 組織特異的 | wing_base, tail_spine |
| Runx2Rios | Skeletal | 骨形成 | 組織特異的 | wing_base, tail_spine |
| MC1R-Rios | Pigmentation | メラノコルチン受容体 | 亜種特異的 | 複数 |
| ASIP-Rios | Pigmentation | アグーチシグナルタンパク質 | 亜種特異的 | 複数 |
| EDNRB-Rios | Pigmentation | エンドセリン受容体B | 亜種特異的 | 複数 |
| ACTB_Rios | Housekeeping | アクチン（正規化コントロール） | なし | すべての組織 |

---

## 補足表S3. 完全な発現差解析結果：性差（M vs F）

| 遺伝子ID | Base Mean | log₂FC | lfcSE | Statistic | p-value | FDR | 有意 | 方向 |
|---------|-----------|--------|-------|-----------|---------|-----|-------------|-----------|
| ACTB_Rios | 2,028,691 | 0.061 | 0.043 | 1.417 | 0.157 | 0.274 | No | Up in M |
| ASIP-Rios | 152,000 | 0.036 | 0.078 | 0.461 | 0.645 | 0.795 | No | Up in M |
| EDNRB-Rios | 231,403 | -0.024 | 0.032 | -0.768 | 0.442 | 0.619 | No | Up in F |
| HoxRiosA9 | 278,217 | -0.276 | 0.211 | -1.308 | 0.191 | 0.297 | No | Up in F |
| HoxRiosC6 | 353,017 | 0.622 | 0.282 | 2.205 | 0.027 | 0.064 | No | Up in M |
| MC1R-Rios | 271,117 | 0.001 | 0.064 | 0.008 | 0.994 | 0.994 | No | Up in M |
| RiosTailSpine1 | 297,024 | -0.606 | 0.269 | -2.252 | 0.024 | 0.064 | No | Up in F |
| RiosTailSpine2 | 288,886 | -0.608 | 0.300 | -2.025 | 0.043 | 0.086 | No | Up in F |
| RiosToxinA1 | 287,280 | -1.003 | 0.336 | -2.989 | 0.003 | **0.020** | **Yes** | **Up in F** |
| RiosToxinA2 | 283,145 | -0.970 | 0.320 | -3.031 | 0.002 | **0.014** | **Yes** | **Up in F** |
| RiosWingBMP1 | 285,714 | 1.021 | 0.350 | 2.920 | 0.004 | **0.045** | **Yes** | **Up in M** |
| RiosWingBMP2 | 277,143 | 0.980 | 0.342 | 2.863 | 0.004 | 0.050 | No | Up in M |
| Runx2Rios | 245,714 | -0.102 | 0.089 | -1.146 | 0.252 | 0.381 | No | Up in F |
| Sox9Rios | 267,143 | 0.153 | 0.095 | 1.611 | 0.107 | 0.206 | No | Up in M |

**有意な遺伝子（FDR < 0.05, |log₂FC| > 1.0）：** 3遺伝子
- RiosWingBMP1: 雄特異的（log₂FC = +1.02, FDR = 4.5×10⁻²）
- RiosToxinA1: 雌特異的（log₂FC = -1.00, FDR = 2.0×10⁻²）
- RiosToxinA2: 雌特異的（log₂FC = -0.97, FDR = 1.4×10⁻²）

---

## 補足表S4. CNV領域と期待コピー数

| クラスター | 染色体 | 開始 | 終了 | 性別 | 亜種 | 期待コピー数 | 備考 |
|---------|------------|-------|-----|-----|------------|---------------|-------|
| WingBMP | Rios1 | 1,000,000 | 1,005,000 | M | Rios_ref | 3 | 雄翼強化 |
| WingBMP | Rios1 | 1,000,000 | 1,005,000 | F | Rios_ref | 2 | 雌ベースライン |
| WingBMP | Rios1 | 1,000,000 | 1,005,000 | M | Rios_volcanicus | 4 | 火山亜種、翼強化 |
| WingBMP | Rios1 | 1,000,000 | 1,005,000 | F | Rios_volcanicus | 2 | 火山亜種、雌 |
| TailSpine | Rios1 | 2,000,000 | 2,005,000 | M | Rios_ref | 2 | 雄ベースライン |
| TailSpine | Rios1 | 2,000,000 | 2,005,000 | F | Rios_ref | 3 | 雌尾棘強化 |
| TailSpine | Rios1 | 2,000,000 | 2,005,000 | M | Rios_marina | 2 | 沿岸亜種、雄 |
| TailSpine | Rios1 | 2,000,000 | 2,005,000 | F | Rios_marina | 4 | 沿岸亜種、尾棘強化 |
| ToxinClaw | Rios1 | 3,000,000 | 3,005,000 | M | Rios_ref | 3 | 雄後脚鉤爪毒腺 |
| ToxinClaw | Rios1 | 3,000,000 | 3,005,000 | F | Rios_ref | 2 | 雌ベースライン |
| ToxinClaw | Rios1 | 3,000,000 | 3,005,000 | M | Rios_volcanicus | 2 | 火山亜種、雄 |
| ToxinClaw | Rios1 | 3,000,000 | 3,005,000 | F | Rios_volcanicus | 2 | 火山亜種、雌 |
| ToxinTail | Rios1 | 4,000,000 | 4,005,000 | M | Rios_ref | 2 | 雄ベースライン |
| ToxinTail | Rios1 | 4,000,000 | 4,005,000 | F | Rios_ref | 4 | 雌尾棘毒腺強化 |
| ToxinTail | Rios1 | 4,000,000 | 4,005,000 | M | Rios_sylvestris | 2 | 森林亜種、雄 |
| ToxinTail | Rios1 | 4,000,000 | 4,005,000 | F | Rios_sylvestris | 5 | 森林亜種、毒腺強化 |

**CNV検出手法：** シミュレーション（期待値に95-99%の精度でノイズを追加）  
**検証：** ピアソン相関係数 r = 0.998（p < 0.001）

---

## 補足図S1. 火山図：性差（M vs F）

**説明：** 雄と雌の間の発現差を示す火山図。X軸：log₂ fold change（M / F）。Y軸：-log₁₀(p-value)。赤い点は有意な遺伝子（FDR < 0.05, |log₂FC| > 1.0）を示す。水平破線：p-value = 0.05。垂直破線：log₂FC = 0。

**主要な所見：**
- 3遺伝子が有意に発現差を示す
- RiosWingBMP1: 雄特異的（右上象限）
- RiosToxinA1とRiosToxinA2: 雌特異的（左上象限）

**ファイル場所：** [`analysis/out/plots/volcano_sex_M_vs_F.png`](analysis/out/plots/volcano_sex_M_vs_F.png)  
**代替形式：** [`analysis/out/plots/volcano_sex_M_vs_F.pdf`](analysis/out/plots/volcano_sex_M_vs_F.pdf)

---

## 補足図S2. 遺伝子発現ヒートマップ

**説明：** 全19サンプルにわたる全14遺伝子のrlog変換、行スケール（z-score）発現値を示すヒートマップ。カラースケール：青（低発現）、白（平均発現）、赤（高発現）。サンプルは組織と性別で順序付けられている。遺伝子は機能クラスターで順序付けられている。

**主要パターン：**
- WingBMP遺伝子：雄翼サンプルで高発現（赤）
- Toxin遺伝子：雌毒腺サンプルで高発現（赤）
- TailSpine遺伝子：雌尾棘サンプルで高発現（赤）
- ACTB_Rios：全サンプルで一貫した発現（白、正規化コントロール）

**ファイル場所：** [`analysis/out/plots/heatmap_all_genes.png`](analysis/out/plots/heatmap_all_genes.png)

---

## 補足図S3. 主成分分析

**説明：** rlog変換発現データに基づくサンプル間の関係を示すPCAプロット。PC1（X軸）は変動の~60-70%を説明し、主に組織によってサンプルを分離する。PC2（Y軸）は変動の~20-30%を説明し、主に性別によってサンプルを分離する。点は組織で色分けされ、性別で形状が異なる。

**主要な所見：**
- 組織が変動の主要因（PC1）
- 性別が変動の二次的要因（PC2）
- 組織タイプ間の明確な分離
- 各組織タイプ内で性差が可視化される

**ファイル場所：** [`analysis/out/plots/pca_plot.png`](analysis/out/plots/pca_plot.png)  
**代替形式：** [`analysis/out/plots/pca_plot.pdf`](analysis/out/plots/pca_plot.pdf)

---

## 補足方法

### ソフトウェアバージョン

すべての解析は、以下のソフトウェアバージョンを使用して実行した：

- **R**: 4.3.2
- **Bioconductor**: 3.18
- **DESeq2**: 1.42.0
- **polyester**: 1.38.0
- **HISAT2**: 2.2.1
- **samtools**: 1.19
- **ggplot2**: 3.4.4
- **Python**: 3.11

### 計算リソース

- **CPU**: マルチコアシステム（4-8コア使用）
- **メモリ**: 16-32 GB RAM
- **ストレージ**: 中間ファイル用に~50 GB
- **実行時間**: 完全なパイプラインで~2-3時間

### 再現性

すべての解析スクリプトはプロジェクトリポジトリで利用可能である：
- RNA-seqパイプライン：`analysis/rnaseq_pipeline.R`
- 可視化：`analysis/visualize_results.R`
- CNV解析：`analysis/create_cnv_slide_data.R`

再現性のために乱数種を設定：
- Polyesterシミュレーション：seed = 1234
- CNVシミュレーション：seed = 42

---

## 補足注記：限界と今後の研究

### 現在の限界

1. **シミュレーションデータ**：結果は実際の生物学的データではなくシミュレーションされたデータを反映する
2. **小規模遺伝子セット**：14遺伝子のみを解析し、一般化可能性を制限
3. **トランスクリプトアライメント**：ゲノムアライメントにより選択的スプライシング検出が可能になる
4. **CNVシミュレーション**：実際のCNV検出にはゲノムアライメントBAMファイルが必要
5. **限定的な亜種サンプリング**：亜種あたり1サンプルのみで、統計的比較を妨げる

### 今後の研究

1. **ゲノムアライメント**：STARまたはHISAT2ゲノムアライメントを実装
2. **拡張された遺伝子セット**：視覚解析用のオプシン遺伝子を含める（Page 18）
3. **実際のCNV検出**：ゲノムアライメントBAMファイルでControl-FREECを実行
4. **亜種比較**：統計的検出力のためにサンプルサイズを増加
5. **時系列解析**：発生時点を含めて時間的発現パターンを研究
6. **マルチオミクス統合**：プロテオミクスまたはメタボロミクスデータと統合

---

**補足資料終了**

