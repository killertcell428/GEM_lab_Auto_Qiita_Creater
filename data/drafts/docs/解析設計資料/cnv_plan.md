# CNV / 亜種メタデータ方針（raw_research_results 準拠）

## 目的
- 性的二形を支えるコピー数差（翼↑♂、尾棘・毒腺↑♀）を in silico で再現し、CNVコーラーで検出可能か検証。
- 亜種（火山=翼強化、森林=毒腺強化、沿岸=尾棘強化）のシナリオを作り、発現・CNVの両面で差を再現。

## ベースコピー数
- `Rios_ref`: WingBMP 2, TailSpine 2, ToxinClaw 2, ToxinTail 2

## 性差の想定
- ♂: WingBMP copy=3, TailSpine=2, ToxinClaw=3（後脚鉤爪の毒腺）, ToxinTail=2
- ♀: WingBMP=2, TailSpine=3, ToxinClaw=2, ToxinTail=4（尾棘の毒腺）
- 注意：最新の生態調査により、オスとメスで毒腺の部位が異なることが判明（オスは後脚鉤爪、メスは尾棘）。そのため、ToxinClawとToxinTailを別クラスターとして扱う。

## 亜種メタ（コピー数想定）
- `Rios_volcanicus` (火山, 翼強化): WingBMP 4, TailSpine 2, ToxinClaw 2, ToxinTail 2
- `Rios_sylvestris` (森林, 毒腺強化): WingBMP 2, TailSpine 2, ToxinClaw 2, ToxinTail 5（メス亜種のため尾棘の毒腺を強化）
- `Rios_marina` (沿岸, 尾棘強化): WingBMP 2, TailSpine 4, ToxinClaw 2, ToxinTail 2

## 実装案
1) CNV BED/VCF を `reference/cnv/` に置き、Control-FREEC で読める形式にする（領域はクラスター単位でまとめる）。
2) RNA-seq シミュレーションではコピー数に比例する期待TPMを設定（expr_matrix_template を亜種・性別で派生）。
3) FASTA/GTF に実座標を追記する場合も、原本は `/home/uecha/data/reference/rios/` に集約し、プロジェクト側はシンボリックリンク。

## Control-FREEC 用メモ
- 入力: BAM, 参照FASTA
- 必要: GC補正, マスキング (optional)
- コンフィグ例: `analysis/freec_config_template.txt`（未作成、テンプレ追加予定）

## 今後のTODO
- CNV領域BEDを作成 (WingBMP, TailSpine, ToxinClaw, ToxinTail)
- 亜種・性差メタデータCSVを作成 (亜種名, sex, copy数, 期待TPM倍率)
- RNA-seqデザインに亜種列を追加し、polyester発現行列を派生させる
- 毒腺の部位特異性を反映した組織サンプル設計（♂後脚鉤爪組織、♀尾棘組織）を追加検討
