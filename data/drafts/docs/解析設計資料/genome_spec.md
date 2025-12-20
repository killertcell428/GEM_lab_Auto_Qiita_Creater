# リオス科ゲノム仕様書（ドラフト）

## 種・亜種・系統（raw_research_results に準拠）
- Genus: Rios
  - R. azureus（リオレウス, ♂翼重視）
  - R. viridis（リオレイア, ♀尾・毒腺重視）
  - 亜種例（役割/形質で区分）
    - R. azureus volcanicus（火山型: 翼強化）
    - R. viridis sylvestris（森林型: 毒腺強化）
    - R. viridis marina（沿岸型: 尾棘強化）

## 性決定様式
- ZW 型（♀: ZW, ♂: ZZ）
- 性差遺伝子ラベル例: dsxRios（dsx相当）、Sox9Rios、FoxL2Rios

## バックボーン
- ベース: Gallus gallus (GRCg7b など最新 Ensembl/NCBI) を `rios_ancestor.fa` として使用予定
- アノテーション: 対応する GTF/GFF3 を `rios_ancestor.gtf` として取得
- 染色体リネーム方針: `chr1`→`Rios1`, `chr2`→`Rios2`, …, `chrZ`→`RiosZ`, `chrW`→`RiosW`, `chrMT`→`RiosMT`

## 既知遺伝子のマッピング（raw_research_results に基づく焦点）
- Hox: HoxA/B/C/D クラスター → HoxRiosA/B/C/D（性差は発現ドメインシフトで再現）
- 軟骨・骨格: Sox9Rios, Runx2Rios（尾椎・翼基部で強調）
- 翼形成: BMPRios2/4/7, WntRios3a/5a/7b（♂で増幅）
- 尾棘形成: HoxRiosA/D 後方ドメイン + Sox9/Runx2 強化（♀で増幅）
- 毒腺: 
  - RiosToxinClawクラスター（♂後脚鉤爪の毒腺、CNVを想定）
  - RiosToxinTailクラスター（♀尾棘の毒腺、CNVを想定）
  - 注意：最新の生態調査により、オスとメスで毒腺の部位が異なることが判明
- 色素: MC1R-Rios, ASIP-Rios, EDNRB-Rios（亜種模様差検証用、♂の翼の炎模様にも関連）
- 視覚・色覚: オプシン遺伝子群（♂で高発現、動体視力・色覚に関与、閃光への脆弱性も関連）
- 哺育器官: RiosJawSpineクラスター（♀下顎の棘、中空構造で哺育器官として機能）

## リオス科特異的挿入（raw_research_results を反映した設計）
- 翼基部軟骨形成（♂主導）: RiosWingBMP1–5（CNV候補）
- 尾椎棘形成（♀主導）: RiosTailSpine1–3（CNV候補）
- 毒腺（部位特異的）:
  - RiosToxinClaw1–3（♂後脚鉤爪の毒腺、CNV前提）
  - RiosToxinTail1–4（♀尾棘の毒腺、CNV前提）
  - 注意：最新の生態調査により、オスとメスで毒腺の部位が異なることが判明したため、別クラスターとして設計
- 哺育器官（♀特異）: RiosJawSpine1–2（♀下顎の棘、中空構造で哺育器官として機能）

## CNV メタデータの想定（性差 + 亜種差）
- 性差:
  - ♂: RiosWingBMP copy↑
  - ♀: RiosTailSpine copy↑, RiosToxin copy↑
- 亜種:
  - volcanicus: 翼クラスター copy↑
  - sylvestris: 毒腺クラスター copy↑
  - marina: 尾棘クラスター copy↑

## ファイル命名
- reference/rios_ancestor.fa, reference/rios_ancestor.gtf
- reference/chrom_rename.tsv（旧→新染色体名）
- reference/rios_genome.fa, reference/rios.gtf（リネーム後）

## 今後の作業
- Ensembl/NCBI からバックボーン取得
- GTFに「Note」「gene_biotype」を付与し、仮想遺伝子を追記
- CNV設定と亜種メタデータを docs/cnv_plan.md に統合
