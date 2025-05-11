#!/bin/bash

# 使用するスレッド数
THREADS=16

# SRA IDリスト
SRA_IDS=("SRR12235286" "SRR12235292" "SRR12235283" "SRR12235284" "SRR12235311" "SRR12235305" "SRR12235314" "SRR12235276" "SRR12235277" "SRR12235291" "SRR12235300" "SRR12235289")

# ダウンロード先ディレクトリ
FASTQ_DIR=~/NGS/20250330Adenoma/fastq

# ログディレクトリ
LOG_DIR=~/NGS/20250330Adenoma/log

# ディレクトリ作成
mkdir -p $FASTQ_DIR $LOG_DIR

# 各SRA IDについてFASTQファイルのダウンロードとgzip圧縮を実行
for SRA_ID in "${SRA_IDS[@]}"; do
    echo "Downloading ${SRA_ID}..."
    
    # ダウンロード
    fasterq-dump --split-files --threads $THREADS $SRA_ID --outdir $FASTQ_DIR
    
    # FASTQファイルをgzip圧縮
    echo "Compressing FASTQ files for ${SRA_ID}..."
    gzip $FASTQ_DIR/${SRA_ID}_1.fastq
    gzip $FASTQ_DIR/${SRA_ID}_2.fastq
    
    echo "Download and compression completed for ${SRA_ID}."
done

echo "All downloads and compressions are completed."

#!/bin/bash

# 使用するスレッド数
THREADS=16

# SRA IDリスト
SRA_IDS=("SRR12235286" "SRR12235292" "SRR12235283" "SRR12235284" "SRR12235311" "SRR12235305" "SRR12235314" "SRR12235276" "SRR12235277" "SRR12235291" "SRR12235300" "SRR12235289")

# 基本ディレクトリ設定
WORK_DIR=~/NGS/20250330Adenoma
FASTQ_DIR=$WORK_DIR/fastq
BAM_DIR=$WORK_DIR/bam
RSEM_DIR=$WORK_DIR/rsem
BIGWIG_DIR=$WORK_DIR/bigwig
LOG_DIR=$WORK_DIR/logs
STAR_REF=~/NGS/ref/STAR_reference   # STARの参照ゲノム
RSEM_REF=~/NGS/ref/RSEM_reference   # RSEMの参照インデックス

# 必要なディレクトリ作成
mkdir -p $FASTQ_DIR $BAM_DIR $RSEM_DIR $BIGWIG_DIR $LOG_DIR

# 各SRA IDに対して解析を実行
for SRA_ID in "${SRA_IDS[@]}"; do
    echo "Starting analysis for ${SRA_ID}..."

    # 各SRA IDに対して解析用ディレクトリを作成
    SAMPLE_BAM_DIR=$BAM_DIR/$SRA_ID
    SAMPLE_RSEM_DIR=$RSEM_DIR/$SRA_ID
    SAMPLE_BIGWIG_DIR=$BIGWIG_DIR/$SRA_ID
    mkdir -p $SAMPLE_BAM_DIR $SAMPLE_RSEM_DIR $SAMPLE_BIGWIG_DIR

    # Step 1: STARアライメント & BAM生成
    echo "Aligning ${SRA_ID} with STAR..."
    ~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR \
        --runThreadN $THREADS \
        --outFilterMultimapNmax 50 \
        --runMode alignReads \
        --genomeDir $STAR_REF \
        --readFilesIn $FASTQ_DIR/${SRA_ID}_1.fastq.gz $FASTQ_DIR/${SRA_ID}_2.fastq.gz \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM \
        --outFileNamePrefix $SAMPLE_BAM_DIR/${SRA_ID}_

    # BAMファイル名
    SORTED_BAM_FILE="${SAMPLE_BAM_DIR}/${SRA_ID}_Aligned.sortedByCoord.out.bam"
    TRANSCRIPTOME_BAM="${SAMPLE_BAM_DIR}/${SRA_ID}_Aligned.toTranscriptome.out.bam"

    # Step 2: BAMのインデックス作成
    echo "Indexing BAM: $SORTED_BAM_FILE"
    samtools index $SORTED_BAM_FILE

    # Step 3: RSEMによるRead Count定量
    echo "Quantifying expression for ${SRA_ID} with RSEM..."
    ~/RSEM-1.3.3/rsem-calculate-expression \
        -p $THREADS \
        --paired-end \
        --bam \
        $TRANSCRIPTOME_BAM \
        $RSEM_REF/RSEM_reference \
        $SAMPLE_RSEM_DIR/${SRA_ID}

    # Step 4: BigWig変換
    echo "Converting BAM to BigWig for ${SRA_ID}..."
    bamCoverage -b $SORTED_BAM_FILE -o $SAMPLE_BIGWIG_DIR/${SRA_ID}.bw --binSize 1 --normalizeUsing CPM -p $THREADS

    echo "Analysis completed for ${SRA_ID}."
done

echo "All analyses completed."
