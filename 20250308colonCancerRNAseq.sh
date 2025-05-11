#!/bin/bash

# 解析用ディレクトリ
WORK_DIR=~/NGS/RNAseq_analysis
FASTQ_DIR=$WORK_DIR/fastq
BAM_DIR=$WORK_DIR/bam
RSEM_DIR=$WORK_DIR/rsem
BIGWIG_DIR=$WORK_DIR/bigwig
LOG_DIR=$WORK_DIR/logs
STAR_REF=~/NGS/ref/STAR_reference   # STARの参照ゲノム
RSEM_REF=~/NGS/ref/RSEM_reference   # RSEMの参照インデックス

# 使用するスレッド数
THREADS=16

# SRA IDと条件名の対応
declare -A SRA_CONDITION_MAP=(
    ["SRR13307407"]="NCM460"
    ["SRR13307410"]="DLD1"
    ["SRR13307411"]="HCT116"
    ["SRR13307412"]="HCT8"
    ["SRR24283812"]="NCM-shGFPrep1"
    ["SRR24283811"]="NCM-shGFPrep2"
    ["SRR24283809"]="NCM-shGFPrep3"
    ["SRR24283780"]="PKO-shGFPrep1"
    ["SRR24283779"]="PKO-shGFPrep2"
    ["SRR24283778"]="PKO-shGFPrep3"
)

# ディレクトリ作成
mkdir -p $FASTQ_DIR $BAM_DIR $RSEM_DIR $BIGWIG_DIR $LOG_DIR

cd $FASTQ_DIR

# Step 1: FASTQのダウンロード
for SRA_ID in "${!SRA_CONDITION_MAP[@]}"; do
    CONDITION=${SRA_CONDITION_MAP[$SRA_ID]}
    PREFIX="${CONDITION}_${SRA_ID}"  # 例: NCM460_SRR13307407

    echo "Downloading ${SRA_ID} (${CONDITION})..."
    
    # ダウンロード（gzipなし）
    fasterq-dump --split-files --threads $THREADS $SRA_ID -O $FASTQ_DIR

    # 圧縮
    echo "Compressing ${PREFIX} FASTQ files..."
    gzip $FASTQ_DIR/${SRA_ID}_1.fastq
    gzip $FASTQ_DIR/${SRA_ID}_2.fastq
done


# Step 2: STARアライメント & BAM生成
for SRA_ID in "${!SRA_CONDITION_MAP[@]}"; do
    CONDITION=${SRA_CONDITION_MAP[$SRA_ID]}
    PREFIX="${CONDITION}_${SRA_ID}"  # 例: NCM460_SRR13307407

    echo "Aligning ${PREFIX} with STAR..."
    
    ~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR \
        --runThreadN $THREADS \
        --outFilterMultimapNmax 50 \
        --runMode alignReads \
        --genomeDir $STAR_REF \
        --readFilesIn $FASTQ_DIR/${SRA_ID}_1.fastq.gz $FASTQ_DIR/${SRA_ID}_2.fastq.gz \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode TranscriptomeSAM \
        --outFileNamePrefix $BAM_DIR/${PREFIX}_

    # BAMファイル名
    SORTED_BAM_FILE="${BAM_DIR}/${PREFIX}_Aligned.sortedByCoord.out.bam"
    TRANSCRIPTOME_BAM="${BAM_DIR}/${PREFIX}_Aligned.toTranscriptome.out.bam"

    # Step 3: BAMのインデックス作成
    echo "Indexing BAM: $SORTED_BAM_FILE"
    samtools index $SORTED_BAM_FILE

    # Step 4: RSEMによるRead Count定量
    echo "Running RSEM for ${PREFIX}..."
    ~/RSEM-1.3.3/rsem-calculate-expression \
        -p $THREADS \
        --paired-end \
        --bam \
        $TRANSCRIPTOME_BAM \
        $RSEM_REF \
        $RSEM_DIR/${PREFIX}

    # Step 5: BigWig変換
    echo "Converting BAM to BigWig for ${PREFIX}..."
    bamCoverage -b $SORTED_BAM_FILE -o $BIGWIG_DIR/${PREFIX}.bw --binSize 1 --normalizeUsing CPM -p $THREADS
done

echo "Pipeline completed successfully!"



##############################################################################################################
#!/bin/bash

# 解析用ディレクトリ
WORK_DIR=~/NGS/RNAseq_analysis
FASTQ_DIR=$WORK_DIR/fastq
BAM_DIR=$WORK_DIR/bam
RSEM_DIR=$WORK_DIR/rsem
BIGWIG_DIR=$WORK_DIR/bigwig
LOG_DIR=$WORK_DIR/logs
STAR_REF=~/NGS/ref/STAR_reference   # STARの参照ゲノム
RSEM_REF=~/NGS/ref/RSEM_reference   # RSEMの参照インデックス

# 使用するスレッド数
THREADS=16

# SRA ID
SRA_ID="SRR24283809"  

# ディレクトリ作成
mkdir -p $FASTQ_DIR $BAM_DIR/$SRA_ID $RSEM_DIR/$SRA_ID $BIGWIG_DIR $LOG_DIR
mkdir -p $RSEM_DIR/${SRA_ID}

# Step 1: FASTQのダウンロード（gzipなし）
echo "Downloading ${SRA_ID}..."
fasterq-dump --split-files --threads $THREADS $SRA_ID

# Step 2: 手動でgzip圧縮
echo "Compressing FASTQ files..."
gzip ${SRA_ID}_1.fastq
gzip ${SRA_ID}_2.fastq

echo "Download and compression completed for ${SRA_ID}."

# Step 2: STARアライメント & BAM生成
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
    --outFileNamePrefix $BAM_DIR/${SRA_ID}/

# BAMファイル名
SORTED_BAM_FILE="${BAM_DIR}/${SRA_ID}/Aligned.sortedByCoord.out.bam"
TRANSCRIPTOME_BAM="${BAM_DIR}/${SRA_ID}/Aligned.toTranscriptome.out.bam"

# Step 3: BAMのインデックス作成
echo "Indexing BAM: $SORTED_BAM_FILE"
samtools index $SORTED_BAM_FILE

# Step 4: RSEMによるRead Count定量
~/RSEM-1.3.3/rsem-calculate-expression \
    -p $THREADS \
    --paired-end \
    --bam \
    $TRANSCRIPTOME_BAM \
    $RSEM_REF/RSEM_reference \
    $RSEM_DIR/${SRA_ID}/${SRA_ID}

# Step 5: BigWig変換
echo "Converting BAM to BigWig for ${SRA_ID}..."
bamCoverage -b $SORTED_BAM_FILE -o $BIGWIG_DIR/${SRA_ID}.bw --binSize 1 --normalizeUsing CPM -p $THREADS

echo "Pipeline completed successfully for ${SRA_ID}!"

#########################################################################################



# RSEMインデックスの作成
~/RSEM-1.3.3/rsem-prepare-reference \
    --gtf ~/NGS/ref/Homo_sapiens.GRCh38.99.gtf \
    ~/NGS/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    ~/NGS/ref/RSEM_reference
# STARインデックスの作成
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR \
    --runMode genomeGenerate \
    --genomeDir ~/NGS/ref/STAR_reference \
    --genomeFastaFiles ~/NGS/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --sjdbGTFfile ~/NGS/ref/Homo_sapiens.GRCh38.99.gtf \
    --runThreadN 16

##################################################################################################
#!/bin/bash

# 解析用ディレクトリ
WORK_DIR=~/NGS/RNAseq_analysis
FASTQ_DIR=$WORK_DIR/fastq
BAM_DIR=$WORK_DIR/bam
RSEM_DIR=$WORK_DIR/rsem
BIGWIG_DIR=$WORK_DIR/bigwig
LOG_DIR=$WORK_DIR/logs
STAR_REF=~/NGS/ref/STAR_reference   # STARの参照ゲノム
RSEM_REF=~/NGS/ref/RSEM_reference   # RSEMの参照インデックス

# 使用するスレッド数
THREADS=16


# SRA ID
SRA_ID="SRR24283809" 

# Step 2: STARアライメント & BAM生成
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
    --outFileNamePrefix $BAM_DIR/${SRA_ID}/

# BAMファイル名
SORTED_BAM_FILE="${BAM_DIR}/${SRA_ID}/Aligned.sortedByCoord.out.bam"
TRANSCRIPTOME_BAM="${BAM_DIR}/${SRA_ID}/Aligned.toTranscriptome.out.bam"

# Step 3: BAMのインデックス作成
echo "Indexing BAM: $SORTED_BAM_FILE"
samtools index $SORTED_BAM_FILE

# Step 4: RSEMによるRead Count定量
~/RSEM-1.3.3/rsem-calculate-expression \
    -p $THREADS \
    --paired-end \
    --bam \
    $TRANSCRIPTOME_BAM \
    $RSEM_REF/RSEM_reference \
    $RSEM_DIR/${SRA_ID}/${SRA_ID}

# Step 5: BigWig変換
echo "Converting BAM to BigWig for ${SRA_ID}..."
bamCoverage -b $SORTED_BAM_FILE -o $BIGWIG_DIR/${SRA_ID}.bw --binSize 1 --normalizeUsing CPM -p $THREADS

echo "Pipeline completed successfully for ${SRA_ID}!"


"SRR13307407"
"SRR13307410"
"SRR13307411"
"SRR13307412"
"SRR24283812"
"SRR24283811"
"SRR24283780"
"SRR24283779"
"SRR24283778"

#!/bin/bash

# 解析用ディレクトリ
WORK_DIR=~/NGS/RNAseq_analysis
FASTQ_DIR=$WORK_DIR/fastq
BAM_DIR=$WORK_DIR/bam
RSEM_DIR=$WORK_DIR/rsem
BIGWIG_DIR=$WORK_DIR/bigwig
LOG_DIR=$WORK_DIR/logs
STAR_REF=~/NGS/ref/STAR_reference   # STARの参照ゲノム
RSEM_REF=~/NGS/ref/RSEM_reference   # RSEMの参照インデックス

# 使用するスレッド数
THREADS=16

# SRA IDのリスト
SRA_IDS=("SRR13307408" "SRR13307409")

# 各SRA IDに対して解析を実行
for SRA_ID in "${SRA_IDS[@]}"; do
    echo "Starting analysis for ${SRA_ID}..."

    # ディレクトリ作成
    mkdir -p $FASTQ_DIR $BAM_DIR/$SRA_ID $RSEM_DIR/$SRA_ID $BIGWIG_DIR $LOG_DIR
    mkdir -p $RSEM_DIR/${SRA_ID}

    # Step 1: FASTQのダウンロード（gzipなし）
    echo "Downloading ${SRA_ID}..."
    fasterq-dump --split-files --threads $THREADS $SRA_ID

    # Step 2: 手動でgzip圧縮
    echo "Compressing FASTQ files..."
    gzip $FASTQ_DIR/${SRA_ID}_1.fastq
    gzip $FASTQ_DIR/${SRA_ID}_2.fastq

    echo "Download and compression completed for ${SRA_ID}."

    # Step 2: STARアライメント & BAM生成
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
        --outFileNamePrefix $BAM_DIR/${SRA_ID}/

    # BAMファイル名
    SORTED_BAM_FILE="${BAM_DIR}/${SRA_ID}/Aligned.sortedByCoord.out.bam"
    TRANSCRIPTOME_BAM="${BAM_DIR}/${SRA_ID}/Aligned.toTranscriptome.out.bam"

    # Step 3: BAMのインデックス作成
    echo "Indexing BAM: $SORTED_BAM_FILE"
    samtools index $SORTED_BAM_FILE

    # Step 4: RSEMによるRead Count定量
    ~/RSEM-1.3.3/rsem-calculate-expression \
        -p $THREADS \
        --paired-end \
        --bam \
        $TRANSCRIPTOME_BAM \
        $RSEM_REF/RSEM_reference \
        $RSEM_DIR/${SRA_ID}/${SRA_ID}

    # Step 5: BigWig変換
    echo "Converting BAM to BigWig for ${SRA_ID}..."
    bamCoverage -b $SORTED_BAM_FILE -o $BIGWIG_DIR/${SRA_ID}.bw --binSize 1 --normalizeUsing CPM -p $THREADS

    echo "Pipeline completed successfully for ${SRA_ID}!"
done

SRR12291410
SRR12291411
SRR12291412
SRR12291413
SRR12291414
SRR12291440
SRR12291442
SRR12291443
SRR12291444
SRR12291445
SRR12291438
SRR12291415
lftp -u 1 <t246013f@yokohama-cu.ac.jp>,<45115414kAisei> ftp://ftp.ncbi.nlm.nih.gov
mput ~/NGS/GEOdeposition/PRO-seq/*
md5sum ~/NGS/GEOdeposition/ChIP-seq/* > ~/NGS/GEOdeposition/ChIP-seq/md5_checksums.txt
md5sum ~/NGS/GEOdeposition/RNA-seq/* > ~/NGS/GEOdeposition/RNA-seq/md5_checksums.txt
md5sum ~/NGS/GEOdeposition/PRO-seq/* > ~/NGS/GEOdeposition/PRO-seq/md5_checksums.txt

put ~/NGS/GEOdeposition/RNA-seq/md5_checksums.txt
put ~/NGS/GEOdeposition/PRO-seq/md5_checksums.txt
put ~/NGS/GEOdeposition/ChIP-seq/md5_checksums.txt
