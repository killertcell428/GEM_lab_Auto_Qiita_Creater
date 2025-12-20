   # 1) バックボーン取得（共通側）
   cd /home/uecha/data/reference/rios
   bash download_backbone.sh

   # 2) 染色体リネーム（共通側で実施、必要なら）
   cd /home/uecha/data/reference/rios
   python3 /home/uecha/work/MonsterGenome_Rios/scripts/rename_chrom.py rios_ancestor.fa  chrom_rename.tsv rios_genome.fa
   python3 /home/uecha/work/MonsterGenome_Rios/scripts/rename_chrom.py rios_ancestor.gtf chrom_rename.tsv rios.gtf

   # 3) プロジェクト側からの参照（シンボリックリンク推奨）
   ln -sf /home/uecha/data/reference/rios/rios_genome.fa /home/uecha/work/MonsterGenome_Rios/reference/rios_genome.fa
   ln -sf /home/uecha/data/reference/rios/rios.gtf      /home/uecha/work/MonsterGenome_Rios/reference/rios.gtf

   # 4) トランスクリプトFASTA準備（ファイル名は仮）
cd /home/uecha/data/reference/rios

python3 - <<'PY'
import random

genes = [
    "RiosWingBMP1","RiosWingBMP2",
    "RiosTailSpine1","RiosTailSpine2",
    "RiosToxinA1","RiosToxinA2",
    "HoxRiosA9","HoxRiosC6",
    "Sox9Rios","Runx2Rios",
    "MC1R-Rios","ASIP-Rios","EDNRB-Rios",
    "ACTB_Rios"
]

def rand_seq(n, gc=0.45):
    # 単純にA/T/G/Cを確率で生成（GC≈45%）
    nts = []
    for _ in range(n):
        nts.append(random.choices(
            population=["A","C","G","T"],
            weights=[(1-gc)/2, gc/2, gc/2, (1-gc)/2],
        )[0])
    return "".join(nts)

lens = {
    "RiosWingBMP1":1500,"RiosWingBMP2":1500,
    "RiosTailSpine1":1500,"RiosTailSpine2":1500,
    "RiosToxinA1":1200,"RiosToxinA2":1200,
    "HoxRiosA9":2000,"HoxRiosC6":2000,
    "Sox9Rios":1400,"Runx2Rios":1400,
    "MC1R-Rios":1100,"ASIP-Rios":1100,"EDNRB-Rios":1100,
    "ACTB_Rios":1800
}

with open("transcripts.fa", "w") as f:
    for g in genes:
        seq = rand_seq(lens[g])
        f.write(f">{g}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
PY

   # 5) リードシミュレーション（scratch へ）
   cd /home/uecha/work/MonsterGenome_Rios
   # 軽量テスト（サンプル4本、リード数少なめ）
   Rscript simulation/polyester/simulate_polyester.R \
     --transcripts /home/uecha/data/reference/rios/transcripts.fa \
     --design simulation/design/design_small.csv \
     --expr simulation/polyester/expr_small.csv \
     --out /home/uecha/scratch/MonsterGenome_Rios/polyester_reads_small \
     --reads_per_sample 1e6 \
     --readlen 100

   # 6) STAR/HISAT2 インデックス作成（共通DB）
   # STAR 例:
   mkdir -p /home/uecha/data/databases/rios/star_index
   STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /home/uecha/data/databases/rios/star_index \
     --genomeFastaFiles /home/uecha/data/reference/rios/rios_genome.fa \
     --sjdbGTFfile /home/uecha/data/reference/rios/rios.gtf \
     --sjdbOverhang 149

   # HISAT2 例:
   mkdir -p /home/uecha/data/databases/rios/hisat2_index
   hisat2-build /home/uecha/data/reference/rios/rios_genome.fa /home/uecha/data/databases/rios/hisat2_index/genome

   # 7) アライメント→featureCounts→DESeq2
   # analysis/rnaseq_pipeline.R の system() コメントを外してから:
   Rscript analysis/rnaseq_pipeline.R \
     --design simulation/design/rios_design.csv \
     --fastq_dir /home/uecha/scratch/MonsterGenome_Rios/polyester_reads \
     --genome /home/uecha/data/reference/rios/rios_genome.fa \
     --gtf /home/uecha/data/reference/rios/rios.gtf \
     --aligner star \
     --threads 8 \
     --outdir analysis/out

   # 8) CNV 用の下ごしらえ（雛形例）
   mkdir -p /home/uecha/work/MonsterGenome_Rios/reference/cnv
   # WingBMP/TailSpine/Toxin の BED/VCF をここに置く
   # Control-FREEC の設定テンプレ（パスは共通ディレクトリに合わせる）:
   cat > /home/uecha/work/MonsterGenome_Rios/analysis/freec_config_template.txt <<'EOF'
   [general]
   bedtools=/usr/bin/bedtools
   samtools=/usr/bin/samtools
   outputDir=/home/uecha/scratch/MonsterGenome_Rios/freec_out
   chrLenFile=/home/uecha/data/reference/rios/rios_genome.fa.fai
   ploidy=2
   window=50000

   [sample]
   mateFile=/path/to/sample.bam
   inputFormat=BAM
   mateOrientation=FR

   [control]
   # optional
   #mateFile=

   [BAF]
   # optional
   EOF