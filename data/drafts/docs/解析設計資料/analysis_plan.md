# 解析手順メモ（目的ドリブン）

信頼情報: `docs/raw_research_results.md`

## ゴール
- 性差: 翼（♂↑）/尾・毒腺（♀↑）の発現差を再現し、HoxシフトとCNV差に整合する DE 結果を得る。
- 毒腺の部位特異性: オス後脚鉤爪の毒腺（ToxinClaw）とメス尾棘の毒腺（ToxinTail）を別クラスターとして解析。
- 亜種: 火山=翼↑, 森林=毒腺↑, 沿岸=尾棘↑ のシナリオを CNV/発現で再現。
- CNV: Control-FREEC などで WingBMP/TailSpine/ToxinClaw/ToxinTail クラスターのコピー差を検出できるか確認。
- 視覚・色覚: オスの優れた動体視力・色覚に関わる遺伝子（オプシン遺伝子群）の発現解析（網膜・視覚野組織）。

## データ配置ポリシー
- 原本リファレンス・インデックス: `/home/uecha/data/reference/rios/`, `/home/uecha/data/databases/rios/`
- プロジェクト側: `/home/uecha/work/MonsterGenome_Rios/reference/` には symlink またはラッパースクリプトのみ。
- 一時展開・中間 BAM/FASTQ: `/home/uecha/scratch/MonsterGenome_Rios/`

## 手順（推奨コマンドの流れ）
1) バックボーン取得（共通側）
   - `bash /home/uecha/data/reference/rios/download_backbone.sh`
   - （必要なら）`python /home/uecha/work/MonsterGenome_Rios/scripts/rename_chrom.py rios_ancestor.fa chrom_rename.tsv rios_genome.fa` を共通側で実行し、プロジェクトは symlink。
2) トランスクリプト FASTA を用意（仮想遺伝子含む）。`simulation/polyester/expr_matrix_template.csv` と ID を合わせる。
3) シミュレーション（polyester）
   - `Rscript simulation/polyester/simulate_polyester.R --transcripts ... --design simulation/design/rios_design.csv --expr simulation/polyester/expr_matrix_template.csv --out /home/uecha/scratch/MonsterGenome_Rios/polyester_reads`
4) アライメント & カウント
   - STAR/HISAT2 インデックスは `/home/uecha/data/databases/rios/` に作成し、`analysis/rnaseq_pipeline.R` のパスを合わせる。
   - `analysis/rnaseq_pipeline.R` の system() コメントを外し、`--design simulation/design/rios_design.csv --fastq_dir /home/uecha/scratch/MonsterGenome_Rios/polyester_reads --outdir analysis/out` などで実行。
5) DE 解析
   - コントラスト例: 
     - M wing_base vs F wing_base（翼性差）
     - F tail_spine vs M tail_spine（尾性差）
     - M claw_toxin vs F claw_toxin（後脚鉤爪の毒腺性差、オス特異的）
     - F tail_toxin vs M tail_toxin（尾棘の毒腺性差、メス特異的）
     - M retina vs F retina（視覚遺伝子の性差、オスで高発現想定）
   - 亜種差: volcanicus vs ref（翼遺伝子群）、sylvestris vs ref（尾棘毒腺群）、marina vs ref（尾棘群）。
6) CNV 解析
   - `reference/cnv/` に WingBMP/TailSpine/ToxinClaw/ToxinTail クラスターの BED/VCF を用意。
   - 注意：最新の生態調査により、オスとメスで毒腺の部位が異なるため、ToxinClaw（後脚鉤爪）とToxinTail（尾棘）を別クラスターとして扱う。
   - Control-FREEC コンフィグを `analysis/freec_config_template.txt` として作成し、BAMとFASTAを指定。出力は `/home/uecha/scratch/MonsterGenome_Rios/freec_out` 推奨。

## 追加TODO
- `reference/cnv/` ディレクトリと BED/VCF の雛形を作成（WingBMP, TailSpine, ToxinClaw, ToxinTail）。
- `analysis/freec_config_template.txt` を追加し、パスを data/databases/rios/ に合わせる。
- STAR/HISAT2 インデックス生成スクリプトを shared/scripts に共通化。
- 毒腺の部位特異性を反映した組織サンプル設計を追加検討（♂後脚鉤爪組織、♀尾棘組織）。
- 視覚・色覚遺伝子の解析のため、網膜・視覚野組織のサンプル設計を追加検討。

