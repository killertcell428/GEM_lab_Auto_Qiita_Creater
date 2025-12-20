# 実験デザイン（raw_research_results に基づく目的ドリブン版）

- 目的: 性的二形の分子基盤（Hox発現シフト + CNV差）を、翼・尾・毒腺で再現する。
- サンプル表: `simulation/design/rios_design.csv`
- 組織: wing_base（翼基部）, tail_spine（尾棘）, toxin_gland（毒腺）, pectoral_muscle, brain, gonad
- 性別: M, F
- リード想定: paired-end 150 bp, 3–5M reads/sample（polyesterで調整）
- 解析ポイント:
  - Sex × Tissue 交互作用 (DESeq2) で Hox/BMP/Wnt（翼） vs HoxA/D/Sox9/Runx2（尾） のシフトを検出
  - 毒腺: F 特異的に Toxin クラスター発現↑（CNVを想定）
  - 亜種拡張: volcanicus=翼↑, sylvestris=毒腺↑, marina=尾棘↑ の発現/コピー差を再現

## 次のステップ
1. バックボーンFASTA/GTF取得（共通側 data/reference/rios/ を原本）
2. 仮想遺伝子座・CNVメタを GTF/BED/メタCSVに反映（cnv_plan に従う）
3. 発現マトリクスを目的に合わせて再スケーリングし polyester でリード生成
4. STAR/HISAT2→featureCounts→DESeq2 で性差・組織差・亜種差が期待どおりか検証
