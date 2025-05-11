
#fastpによるtriming 
FASTQ_DIR="./"
OUTPUT_DIR="./fastp_output"
LOG_DIR="./fastp_logs"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

for FILE in ${FASTQ_DIR}/*.fq.gz ${FASTQ_DIR}/*.fastq.gz; do
  BASENAME=$(basename "${FILE}" | sed 's/_1\.fq\.gz//; s/_2\.fq\.gz//; s/\.fastq\.gz//')

  # paired-endの場合
  if [[ -f "${FASTQ_DIR}/${BASENAME}_1.fq.gz" && -f "${FASTQ_DIR}/${BASENAME}_2.fq.gz" ]]; then
    fastp \
      -i "${FASTQ_DIR}/${BASENAME}_1.fq.gz" \
      -I "${FASTQ_DIR}/${BASENAME}_2.fq.gz" \
      -o "${OUTPUT_DIR}/${BASENAME}_1_trimmed.fq.gz" \
      -O "${OUTPUT_DIR}/${BASENAME}_2_trimmed.fq.gz" \
      --html "${LOG_DIR}/${BASENAME}_fastp.html" \
      --json "${LOG_DIR}/${BASENAME}_fastp.json" \
      --thread 16
  # single-endの場合
  elif [[ -f "${FASTQ_DIR}/${BASENAME}.fastq.gz" ]]; then
    fastp \
      -i "${FASTQ_DIR}/${BASENAME}.fastq.gz" \
      -o "${OUTPUT_DIR}/${BASENAME}_trimmed.fq.gz" \
      --html "${LOG_DIR}/${BASENAME}_fastp.html" \
      --json "${LOG_DIR}/${BASENAME}_fastp.json" \
      --thread 16
  else
    echo "Skipping ${FILE}: Matching files not found for paired-end or single-end processing."
  fi
done

echo "fastp trimming completed!"

#samファイルの作成　マッピング

# ディレクトリ設定
FASTQ_DIR="./fastp_output" # トリミング済みFastQファイルのディレクトリ
OUTPUT_DIR="./bowtie2_output" # Bowtie2マッピング結果の保存ディレクトリ
INDEX_DIR="/home/shintaro/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"

# 出力ディレクトリを作成
mkdir -p "${OUTPUT_DIR}"

# FastQファイルの処理
for FILE in ${FASTQ_DIR}/*_trimmed.fq.gz; do
  BASENAME=$(basename "${FILE}" | sed 's/_1_trimmed\.fq\.gz//; s/_2_trimmed\.fq\.gz//; s/_trimmed\.fq\.gz//')
  
  # ペアエンドの場合
  if [[ -f "${FASTQ_DIR}/${BASENAME}_1_trimmed.fq.gz" && -f "${FASTQ_DIR}/${BASENAME}_2_trimmed.fq.gz" ]]; then
    bowtie2 \
      -x ${INDEX_DIR} \
      -1 "${FASTQ_DIR}/${BASENAME}_1_trimmed.fq.gz" \
      -2 "${FASTQ_DIR}/${BASENAME}_2_trimmed.fq.gz" \
      -S "${OUTPUT_DIR}/${BASENAME}.sam" \
      --very-sensitive \
      --threads 16
  # シングルエンドの場合
  elif [[ -f "${FASTQ_DIR}/${BASENAME}_trimmed.fq.gz" ]]; then
    bowtie2 \
      -x ${INDEX_DIR} \
      -U "${FASTQ_DIR}/${BASENAME}_trimmed.fq.gz" \
      -S "${OUTPUT_DIR}/${BASENAME}.sam" \
      --very-sensitive \
      --threads 16
  else
    echo "Skipping ${FILE}: Matching files not found for paired-end or single-end processing."
  fi
done

echo "Bowtie2 mapping completed!"

#bamファイルの作成-q 15バージョン

# スレッド数
THREADS=16

# 入力と出力ディレクトリ
SAM_DIR="./bowtie2_output"
BAM_DIR="./bam_output"

# BAM出力ディレクトリを作成（存在しない場合）
mkdir -p ${BAM_DIR}

# 処理を実行
for FILE in ${SAM_DIR}/*.sam; do
    # サンプル名を取得
    BASENAME=$(basename ${FILE} .sam)
    
    # BAMファイル名
    BAM_FILE=${BAM_DIR}/${BASENAME}.bam
    
    # BAMファイルに変換（-q 15 フィルタ適用）
    samtools view -@ ${THREADS} -bhS -F 0x4 -q 15 ${FILE} | \
    samtools sort -@ ${THREADS} -o ${BAM_FILE}
    
    # インデックス作成
    samtools index -@ ${THREADS} ${BAM_FILE}
    
    echo "Processed ${BAM_FILE}"
done

#bamファイルの作成-q 0バージョン

# スレッド数
THREADS=16

# 入力と出力ディレクトリ
SAM_DIR="./bowtie2_output"
BAM_DIR="./bam_output_1"

# BAM出力ディレクトリを作成（存在しない場合）
mkdir -p ${BAM_DIR}

# 処理を実行
for FILE in ${SAM_DIR}/*.sam; do
    # サンプル名を取得
    BASENAME=$(basename ${FILE} .sam)
    
    # BAMファイル名
    BAM_FILE=${BAM_DIR}/${BASENAME}.bam
    
    # BAMファイルに変換（-q 0 フィルタ適用）
    samtools view -@ ${THREADS} -bhS -F 0x4 -q 0 ${FILE} | \
    samtools sort -@ ${THREADS} -o ${BAM_FILE}
    
    # インデックス作成
    samtools index -@ ${THREADS} ${BAM_FILE}
    
    echo "Processed ${BAM_FILE}"
done


#bigwigファイルの作成

# 入力と出力ディレクトリ
BAM_DIR="./bam_output_1"
BW_DIR="./bigwig_output"

# BigWig 出力ディレクトリを作成（存在しない場合）
mkdir -p ${BW_DIR}

# BAMファイルからBigWig生成
for BAM_FILE in ${BAM_DIR}/*.bam; do
    # サンプル名取得
    BASENAME=$(basename ${BAM_FILE} .bam)
    
    # BigWig出力ファイル
    BW_FILE=${BW_DIR}/${BASENAME}.bw
    
    # BigWig生成
    bamCoverage -p 16 -bs 1 -b ${BAM_FILE} -o ${BW_FILE} -of bigwig --normalizeUsing CPM
    
    echo "Generated ${BW_FILE}"
done


# Coilinのピークコール
macs2 callpeak -t ~/NGS/20241127_GFP_coilin_ChIP/bam_output_1/SRR1390685.bam -c ~/NGS/20241127_GFP_coilin_ChIP/bam_output_1/SRR1390686.bam -f BAM -g hs -n Coilin_ChIP --outdir ~/NGS/20241127_GFP_coilin_ChIP/macs2_output_1
macs2 callpeak -t ~/NGS/20241127_GFP_coilin_ChIP/bam_output_1/SRR1390685.bam -c ~/NGS/20241127_GFP_coilin_ChIP/bam_output_1/SRR1390686.bam -f BAM -g hs -n Coilin_ChIP --outdir ~/NGS/20241127_GFP_coilin_ChIP/macs2_output_1

# DMSOのピークコール
macs2 callpeak -t ~/NGS/20241127_GFP_coilin_ChIP/bam_output/DMSO_1.bam -t ~/NGS/20241127_GFP_coilin_ChIP/bam_output/DMSO_2.bam -c ~/NGS/20241127_GFP_coilin_ChIP/bam_output/Input.bam -f BAM -g hs -n DMSO_ChIP --outdir ~/NGS/20241127_GFP_coilin_ChIP/macs2_output



awk '$3 == "gene" && $0 ~ /protein_coding/ {split($9, a, ";"); for(i in a) if(a[i] ~ /gene_id/) gene_id=substr(a[i], 9, length(a[i])-8); print $1"\t"$4-1"\t"$5"\t"gene_id"\t0\t"$7}' ~/NGS/gencode.v45.annotation.gtf > hg38_protein_coding_gene_with_strand_and_name.bed
awk '$3 == "gene" && $0 ~ /snRNA/ {
    match($0, /gene_name "([^"]+)"/, name);
    strand = $7;
    print $1 "\t" $4-1 "\t" $5 "\t" name[1] "\t0\t" strand
}' ~/NGS/gencode.v45.annotation.gtf > hg38_snRNA_with_strand_and_name.bed

awk '$3 == "gene" && $0 ~ /snRNA/ {split($9, a, ";"); for(i in a) if(a[i] ~ /gene_id/) gene_id=substr(a[i], 9, length(a[i])-8); print $1"\t"$4-1"\t"$5"\t"gene_id"\t0\t"$7}' ~/NGS/gencode.v45.annotation.gtf > hg38_snRNA_with_strand_and_name.bed
##major snRNAのみ
"\\wsl.localhost\Ubuntu\home\shintaro\NGS\20250115ICE1AFF4_Pol2ChIP\fastq\bigwig_output\ICE1_pol2ChIP2.bw"
cd ~/NGS/20250115ICE1AFF4_Pol2ChIP/fastq/bigwig_output
computeMatrix scale-regions -p 16 \
  -S AFF4_pol2ChIP1.bw AFF4_pol2ChIP2.bw ICE1_pol2ChIP1.bw ICE1_pol2ChIP2.bw \
  -R hg38_majorRNU.bed \
  --outFileName ~/NGS/20250115ICE1AFF4_Pol2ChIP/fastq/bigwig_output/ICE1AFF4.txt.gz \
  -a 1000 -b 1000 --skipZeros

plotProfile \
  --matrixFile ~/NGS/20250115ICE1AFF4_Pol2ChIP/fastq/bigwig_output/ICE1AFF4.txt.gz \
  --outFileName ~/NGS/20250115ICE1AFF4_Pol2ChIP/fastq/bigwig_output/MetaGenePlot.pdf \
  --outFileNameData ~/NGS/20250115ICE1AFF4_Pol2ChIP/fastq/bigwig_output/MetaGenePlot_Data.txt \
  --regionsLabel "snRNA" \
  --samplesLabel "AFF4r1" "AFF4r2" "ICE1r1" "ICE1r2" \
  --perGroup \
  --legendLocation "best" \
  --colors blue red green black \

cd ~/NGS/20250115ICE1AFF4_Pol2ChIP/fastq/bigwig_output
computeMatrix scale-regions -p 16 \
  -S AFF4_pol2ChIP2.bw ICE1_pol2ChIP2.bw \
  -R hg38_majorRNU.bed \
  --outFileName ~/NGS/20250115ICE1AFF4_Pol2ChIP/fastq/bigwig_output/ICE1AFF4_r2.txt.gz \
  -a 1000 -b 1000 --skipZeros

computeMatrix scale-regions -p 16 \
  -S ICE1_pol2ChIP2.bw \
  -R hg38_majorRNU.bed \
  --outFileName ~/NGS/20250115ICE1AFF4_Pol2ChIP/fastq/bigwig_output/ICE1_r2.txt.gz \
  -a 1000 -b 1000 --skipZeros
computeMatrix scale-regions -p 16 \
  -S AFF4_pol2ChIP2.bw \
  -R hg38_majorRNU.bed \
  --outFileName ~/NGS/20250115ICE1AFF4_Pol2ChIP/fastq/bigwig_output/AFF4_r2.txt.gz \
  -a 1000 -b 1000 --skipZeros
gunzip ICE1_r2.txt.gz
gunzip AFF4_r2.txt.gz

plotProfile \
  --matrixFile ~/NGS/20250115ICE1AFF4_Pol2ChIP/fastq/bigwig_output/ICE1AFF4_r2.txt.gz \
  --outFileName ~/NGS/20250115ICE1AFF4_Pol2ChIP/fastq/bigwig_output/MetaGenePlot_r2.pdf \
  --outFileNameData ~/NGS/20250115ICE1AFF4_Pol2ChIP/fastq/bigwig_output/MetaGenePlot_Data_r2.txt \
  --regionsLabel "snRNA" \
  --samplesLabel "AFF4" "ICE1" \
  --perGroup \
  --legendLocation "best" \
    --yMin 0
  
gunzip ICE1AFF4_r2.txt.gz
##snRNAすべて(RNU6系列は除いてある)
computeMatrix scale-regions -p 16 \
  -S DMSO_average.bw Input.bw \
  -R hg38_snRNA_filtered.bed hg38_protein_coding_gene_with_strand_and_name.bed 221120_hist.bed \
  --outFileName ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/ICE2.txt_1.gz \
  -a 1000 -b 1000 --skipZeros
plotProfile \
  --matrixFile ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/ICE2.txt_1.gz \
  --outFileName ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/MetaGenePlot_1.pdf \
  --outFileNameData ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/MetaGenePlot_Data_1.txt \
  --regionsLabel "snRNA" "Histone" "ProteinCoding" \
  --samplesLabel "DMSO_Average" "Input" \
  --perGroup \
  --legendLocation "best" \
  --colors blue red \
##snRNA
computeMatrix scale-regions -p 16 \
  -S DMSO_average.bw Input.bw \
  -R hg38_snRNA_filtered.bed \
  --outFileName ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/ICE2.txt_snRNA.gz \
  -a 1000 -b 1000 --skipZeros
plotProfile \
  --matrixFile ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/ICE2.txt_snRNA.gz \
  --outFileName ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/MetaGenePlot_snRNA.pdf \
  --outFileNameData ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/MetaGenePlot_Data_snRNA.txt \
  --regionsLabel "snRNA" \
  --samplesLabel "DMSO_Average" "Input" \
  --perGroup \
  --legendLocation "best" \
  --colors blue red \

##Protein
computeMatrix scale-regions -p 16 \
  -S DMSO_average.bw Input.bw \
  -R hg38_protein_coding_gene_with_strand_and_name.bed \
  --outFileName ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/ICE2.txt_pro.gz \
  -a 1000 -b 1000 --skipZeros
plotProfile \
  --matrixFile ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/ICE2.txt_pro.gz \
  --outFileName ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/MetaGenePlot_pro.pdf \
  --outFileNameData ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/MetaGenePlot_Data_pro.txt \
  --regionsLabel "ProteinCoding" \
  --samplesLabel "DMSO_Average" "Input" \
  --perGroup \
  --legendLocation "best" \
  --colors blue red \

##histone
computeMatrix scale-regions -p 16 \
  -S DMSO_average.bw Input.bw \
  -R 221120_hist.bed \
  --outFileName ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/ICE2.txt_hist.gz \
  -a 1000 -b 1000 --skipZeros
plotProfile \
  --matrixFile ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/ICE2.txt_hist.gz \
  --outFileName ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/MetaGenePlot_hist.pdf \
  --outFileNameData ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/MetaGenePlot_Data_hist.txt \
  --regionsLabel "Histone" \
  --samplesLabel "DMSO_Average" "Input" \
  --perGroup \
  --legendLocation "best" \
  --colors blue red \

plotHeatmap \
  --matrixFile ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/ICE2.txt_1.gz \
  --outFileName ~/NGS/20241127_GFP_coilin_ChIP/bigwig_output/Heatmap_1.pdf \
  --regionsLabel "snRNA" "ProteinCoding" "Histone" \
  --samplesLabel "DMSO_Average" "Input" \
  --colorMap "RdBu" \
  --heatmapHeight 10 \
  --heatmapWidth 5 \
  --legendLocation "best" \
  --sortRegions descend \
  --zMin -2 \
  --zMax 2