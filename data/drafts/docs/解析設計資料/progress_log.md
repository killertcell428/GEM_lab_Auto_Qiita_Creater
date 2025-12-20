# 解析進捗ログ

## 完了したステップ（2025-12-11）

### 1. バックボーン取得・リネーム ✅
- **場所**: `/home/uecha/data/reference/rios/`
- **取得ファイル**:
  - `rios_ancestor.fa` (Gallus gallus GRCg7b, Ensembl release 112)
  - `rios_ancestor.gtf` (同上)
- **リネーム後**:
  - `rios_genome.fa` (chr1→Rios1, chrZ→RiosZ, chrW→RiosW, chrMT→RiosMT など)
  - `rios.gtf` (同様に染色体名リネーム)
- **プロジェクト側リンク**:
  - `/home/uecha/work/MonsterGenome_Rios/reference/rios_genome.fa` → `/home/uecha/data/reference/rios/rios_genome.fa`
  - `/home/uecha/work/MonsterGenome_Rios/reference/rios.gtf` → `/home/uecha/data/reference/rios/rios.gtf`

### 2. トランスクリプトFASTA作成 ✅
- **場所**: `/home/uecha/data/reference/rios/transcripts.fa`
- **内容**: 14遺伝子の合成配列（RiosWingBMP1/2, RiosTailSpine1/2, RiosToxinA1/2, HoxRiosA9/C6, Sox9Rios, Runx2Rios, MC1R-Rios, ASIP-Rios, EDNRB-Rios, ACTB_Rios）
- **作成方法**: Python3 でランダム配列生成（GC含量≈45%）

### 3. RNA-seqリードシミュレーション（polyester）✅
- **環境**: conda `rios-r43` (R 4.3 + Bioconductor 3.18)
- **入力**:
  - トランスクリプト: `/home/uecha/data/reference/rios/transcripts.fa`
  - デザイン: `simulation/design/rios_design.csv` (19サンプル)
  - 発現行列: `simulation/polyester/expr_matrix_template.csv`
- **出力**: `/home/uecha/scratch/MonsterGenome_Rios/polyester_reads/`
- **パラメータ**: ペアエンド 150bp, サンプルあたり 5M reads（デフォルト）

---

## デバッグで分かったこと

### polyester 実行時の問題と解決策

#### 1. R環境のバージョン不一致
- **問題**: `rios-env` (R 4.4) では `polyester` が Bioconductor 3.20 に未対応
- **解決**: R 4.3 専用環境 `rios-r43` を作成
  ```bash
  conda create -y -n rios-r43 -c conda-forge -c bioconda \
    r-base=4.3 bioconductor-polyester=1.38.0 ...
  ```

#### 2. パッケージ依存
- **必須**: `Biostrings` (readDNAStringSet 用), `polyester`, `data.table`, `optparse`
- **注意**: `GenomeInfoDbData` の post-link スクリプトが conda で失敗する場合あり → R から直接インストール可

#### 3. simulate_experiment の引数
- **必須引数**:
  - `fasta`: ファイルパス文字列（DNAStringSet は内部で読み込まれる）
  - `reads_per_transcript`: 行列（行=遺伝子, 列=サンプル）
  - `fold_changes`: 同サイズの行列（全1でも可）
  - `num_reps`: サンプル列数（明示指定推奨）
- **ID整合性**: FASTAヘッダと expr_matrix の gene_id が完全一致が必要

#### 4. メモリ消費
- **問題**: 19サンプル×14遺伝子×5M reads でメモリ消費が大きい
- **対策**: 軽量版（4サンプル, 1M reads, readlen 100bp）を用意
  - `simulation/design/design_small.csv`
  - `simulation/polyester/expr_small.csv`
  - 実行例: `docs/memo.sh` 参照

---

## 次回の続き手順

### 4. STAR/HISAT2 インデックス作成（未着手）
- **場所**: `/home/uecha/data/databases/rios/star_index`, `/home/uecha/data/databases/rios/hisat2_index`
- **STAR コマンド例**:
  ```bash
  mkdir -p /home/uecha/data/databases/rios/star_index
  STAR --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir /home/uecha/data/databases/rios/star_index \
    --genomeFastaFiles /home/uecha/data/reference/rios/rios_genome.fa \
    --sjdbGTFfile /home/uecha/data/reference/rios/rios.gtf \
    --sjdbOverhang 149
  ```
- **HISAT2 コマンド例**:
  ```bash
  mkdir -p /home/uecha/data/databases/rios/hisat2_index
  hisat2-build /home/uecha/data/reference/rios/rios_genome.fa \
    /home/uecha/data/databases/rios/hisat2_index/genome
  ```

### 5. アライメント→featureCounts→DESeq2（未着手）
- **スクリプト**: `analysis/rnaseq_pipeline.R`
- **注意**: `system()` コメントを外し、パスを共通ディレクトリ方針に合わせる
- **入力**: `/home/uecha/scratch/MonsterGenome_Rios/polyester_reads/` の FASTQ
- **出力**: `analysis/out/` (BAM, gene_counts.txt, dds.rds, vst.rds)

### 6. CNV検出（未着手）
- **ツール**: Control-FREEC
- **設定**: `analysis/freec_config_template.txt` を作成済み（要編集）
- **CNV領域定義**: `reference/cnv/` に BED/VCF を配置予定

---

## ファイル構成（現在）

```
/home/uecha/
├── data/
│   ├── reference/rios/
│   │   ├── rios_ancestor.fa, rios_ancestor.gtf (元データ)
│   │   ├── rios_genome.fa, rios.gtf (リネーム後)
│   │   ├── transcripts.fa (14遺伝子の合成配列)
│   │   ├── chrom_rename.tsv (リネーム対応表)
│   │   └── download_backbone.sh (取得スクリプト)
│   └── databases/rios/ (インデックス用、未作成)
│       ├── star_index/
│       └── hisat2_index/
├── work/MonsterGenome_Rios/
│   ├── reference/ (シンボリックリンク)
│   ├── simulation/
│   │   ├── design/
│   │   │   ├── rios_design.csv (19サンプル)
│   │   │   └── design_small.csv (4サンプル軽量版)
│   │   └── polyester/
│   │       ├── simulate_polyester.R (修正済み)
│   │       ├── expr_matrix_template.csv (19サンプル)
│   │       └── expr_small.csv (4サンプル軽量版)
│   └── analysis/
│       └── rnaseq_pipeline.R (system() コメントアウト中)
└── scratch/MonsterGenome_Rios/
    └── polyester_reads/ (FASTQ生成済み)
```

---

## 軽量版の使い方（メモリ節約）

```bash
conda run -n rios-r43 Rscript /home/uecha/work/MonsterGenome_Rios/simulation/polyester/simulate_polyester.R \
  --transcripts /home/uecha/data/reference/rios/transcripts.fa \
  --design /home/uecha/work/MonsterGenome_Rios/simulation/design/design_small.csv \
  --expr /home/uecha/work/MonsterGenome_Rios/simulation/polyester/expr_small.csv \
  --out /home/uecha/scratch/MonsterGenome_Rios/polyester_reads_small \
  --reads_per_sample 1e6 \
  --readlen 100
```

---

## 次回作業開始時のチェックリスト

- [ ] conda 環境 `rios-r43` が有効か確認
- [ ] `/home/uecha/scratch/MonsterGenome_Rios/polyester_reads/` に FASTQ が存在するか確認
- [ ] STAR/HISAT2 が `rios-r43` 環境にインストール済みか確認（なければ追加）
- [ ] `analysis/rnaseq_pipeline.R` の `system()` コメントを外し、パスを確認
- [ ] インデックス作成 → アライメント → DE の順で実行

---

最終更新: 2025-12-11

---

## 完了したステップ（2025-12-13）

### 4. トランスクリプトインデックス作成 ✅
- **場所**: `/home/uecha/data/databases/rios/hisat2_transcript_index/`
- **理由**: polyesterで生成したリードは`transcripts.fa`由来のため、ゲノムインデックスではアライメントできない
- **作成方法**:
  ```bash
  hisat2-build -p 4 \
    /home/uecha/data/reference/rios/transcripts.fa \
    /home/uecha/data/databases/rios/hisat2_transcript_index/transcripts
  ```
- **結果**: インデックスファイル8個（合計約4.2MB）作成完了

### 5. RNA-seqパイプライン実行 ✅
- **スクリプト**: `analysis/rnaseq_pipeline.R`
- **環境**: conda `rios-env` (R 4.3 + DESeq2)
- **実行内容**:
  1. **アライメント**: HISAT2で19サンプルすべてアライメント（99.99%アライメント率）
  2. **カウント**: トランスクリプトIDベースのカウント（`count_transcripts()`関数）
  3. **DESeq2解析**: 完了（デザイン式: `sex + tissue`）
- **出力ファイル**:
  - `analysis/out/gene_counts.txt` (14遺伝子×19サンプル)
  - `analysis/out/dds.rds` (DESeq2オブジェクト)
  - `analysis/out/vst.rds` (rlog変換オブジェクト、遺伝子数が少ないためVSTではなくrlogを使用)

### 6. パイプライン修正・デバッグ ✅

#### 問題1: アライメント率0%
- **原因**: ゲノムインデックスに対してトランスクリプトリードをアライメントしていた
- **解決**: トランスクリプトインデックスを作成し、自動判定機能を追加

#### 問題2: featureCountsでカウント0%
- **原因**: トランスクリプトアライメントのBAMにはゲノム座標がないため、ゲノムGTFとマッチしない
- **解決**: `count_transcripts()`関数を作成し、BAMから直接トランスクリプトIDをカウント

#### 問題3: DESeq2のモデル行列がフルランクでない
- **原因**: 亜種サンプルが各1つしかなく、`subspecies:sex:tissue`の交互作用が計算できない
- **解決**: 亜種サンプル数をチェックし、少ない場合は`sex + tissue`の主効果のみを使用

#### 問題4: VST変換エラー
- **原因**: 遺伝子数が少ない（14遺伝子）ため、`vst()`関数が使用できない
- **解決**: 遺伝子数が1000未満の場合は自動的に`rlog()`を使用

---

## 次回の続き手順

### 7. DESeq2結果の可視化（未着手）
- **火山図**: 性差（M vs F）、組織差（wing_base vs tail_spine など）
- **ヒートマップ**: 14遺伝子の発現パターン
- **PCA**: サンプル間の関係性
- **スクリプト**: `analysis/visualize_results.R` を作成予定

### 8. CNV検出（未着手）
- **ツール**: Control-FREEC
- **設定**: `analysis/freec_config_template.txt` を作成済み（要編集）
- **CNV領域定義**: `reference/cnv/` に BED/VCF を配置予定
- **入力**: BAMファイル（`analysis/out/*.bam`）

---

## ファイル構成（現在）

```
/home/uecha/
├── data/
│   ├── reference/rios/
│   │   ├── rios_ancestor.fa, rios_ancestor.gtf (元データ)
│   │   ├── rios_genome.fa, rios.gtf (リネーム後)
│   │   ├── transcripts.fa (14遺伝子の合成配列)
│   │   ├── chrom_rename.tsv (リネーム対応表)
│   │   └── download_backbone.sh (取得スクリプト)
│   └── databases/rios/
│       ├── star_index/ (未作成、メモリ不足のため保留)
│       ├── hisat2_index/ (ゲノム用、未使用)
│       └── hisat2_transcript_index/ ✅ (トランスクリプト用、作成済み)
├── work/MonsterGenome_Rios/
│   ├── reference/ (シンボリックリンク)
│   ├── simulation/
│   │   ├── design/
│   │   │   └── rios_design.csv (19サンプル)
│   │   └── polyester/
│   │       ├── simulate_polyester.R
│   │       └── expr_matrix_template.csv
│   ├── analysis/
│   │   ├── rnaseq_pipeline.R ✅ (修正完了、動作確認済み)
│   │   └── out/ ✅
│   │       ├── *.bam (19サンプルのアライメント結果)
│   │       ├── gene_counts.txt (カウント結果)
│   │       ├── dds.rds (DESeq2オブジェクト)
│   │       └── vst.rds (rlog変換オブジェクト)
│   └── scripts/
│       ├── rename_chrom.py
│       ├── build_index.sh
│       ├── build_transcript_index.sh ✅ (新規作成)
│       └── convert_fasta_to_fastq.py
└── scratch/MonsterGenome_Rios/
    └── polyester_reads/ (FASTA形式、自動でFASTQに変換)
```

---

## 次回作業開始時のチェックリスト

- [x] トランスクリプトインデックス作成済み
- [x] アライメント完了（19サンプル、99.99%アライメント率）
- [x] カウント完了（14遺伝子）
- [x] DESeq2解析完了
- [ ] DESeq2結果の可視化（火山図、ヒートマップ、PCA）
- [ ] CNV検出（Control-FREEC）

---

最終更新: 2025-12-13

