
# 変数定義
BOWTIE_INDEX=~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index
CHROM_SIZES=hg38.chrom.sizes
THREADS=16  # 並列実行スレッド数

mkdir -p output trimmed_logs bam bw

# サンプルごとに処理
for R1 in *_R1_001.fastq.gz; do
    SAMPLE=$(basename $R1 _R1_001.fastq.gz)
    R2=${SAMPLE}_R2_001.fastq.gz

    echo "Processing: $SAMPLE"

    ## 1. アダプタートリミング & クオリティチェック
    fastp -i $R1 -I $R2 -o trimmed_logs/${SAMPLE}_R1_trimmed.fastq.gz -O trimmed_logs/${SAMPLE}_R2_trimmed.fastq.gz \
          --detect_adapter_for_pe -h trimmed_logs/${SAMPLE}_fastp.html -j trimmed_logs/${SAMPLE}_fastp.json \
          -w $THREADS

    ## 2. アライメント (bowtie2)
    bowtie2 -x $BOWTIE_INDEX -1 trimmed_logs/${SAMPLE}_R1_trimmed.fastq.gz -2 trimmed_logs/${SAMPLE}_R2_trimmed.fastq.gz \
            -S output/${SAMPLE}.sam --threads $THREADS --no-mixed --no-discordant

    ## 3. BAM処理（ソート & 重複除去）
    samtools view -bS output/${SAMPLE}.sam | samtools sort -o bam/${SAMPLE}.sorted.bam
    samtools index bam/${SAMPLE}.sorted.bam
    rm output/${SAMPLE}.sam  # SAM削除してディスク節約

done

echo "Pipeline completed!"


#!/bin/bash

# 変数定義
BOWTIE_INDEX=~/bowtie2_index/BDGP6/BDGP6  # ハエのBowtie2インデックス
CHROM_SIZES=dm6.chrom.sizes  # ハエのchrom.sizes
THREADS=16  # 並列実行スレッド数


# サンプルごとに処理
for R1 in ~/NGS/250406PROseq_re/*_R1_001.fastq.gz; do
    SAMPLE=$(basename $R1 _R1_001.fastq.gz)
    R2=~/NGS/250406PROseq_re/${SAMPLE}_R2_001.fastq.gz

    echo "Processing: $SAMPLE (mapping to dm6)"

    ## 2. アライメント (bowtie2)
    bowtie2 -x $BOWTIE_INDEX \
            -1 trimmed_logs/${SAMPLE}_R1_trimmed.fastq.gz \
            -2 trimmed_logs/${SAMPLE}_R2_trimmed.fastq.gz \
            -S output_dm6/${SAMPLE}.sam \
            --threads $THREADS --no-mixed --no-discordant

    ## 3. BAM処理（ソート & 重複除去）
    samtools view -bS output_dm6/${SAMPLE}.sam | samtools sort -o bam_dm6/${SAMPLE}.sorted.bam
    samtools index bam_dm6/${SAMPLE}.sorted.bam
    rm output_dm6/${SAMPLE}.sam  # SAM削除してディスク節約

done
echo "Pipeline for dm6 mapping completed!"

#!/bin/bash

# 変数定義
CHROM_SIZES=hg38.chrom.sizes
THREADS=16
mkdir -p bw_readLength_normalized

echo "Step 1: Collecting read length distributions from Drosophila BAMs"

# ハエのリード長ごとのリード数を取得
for SAMPLE in bam_dm6/*.sorted.bam; do
    samtools view "$SAMPLE" | awk '{print length($10)}' | sort | uniq -c | awk '{print $2, $1}' > "${SAMPLE%.sorted.bam}_read_length.txt"
done

# すべてのサンプルのリード長ごとのカウントを統合
echo "Step 2: Calculating median read count per length"

paste bam_dm6/*_read_length.txt | awk '
{
    for (i=2; i<=NF; i+=2) {
        read_length[$1][i/2] = $i;
    }
}
END {
    for (len in read_length) {
        n = asort(read_length[len]);
        median = (n % 2 == 1) ? read_length[len][int(n/2)+1] : (read_length[len][n/2] + read_length[len][n/2+1]) / 2;
        print len, median;
    }
}' > dm6_median_read_length.txt

# サンプルごとに正規化を適用
echo "Step 3: Normalizing human reads using Drosophila medians"

for SAMPLE in bam/*.sorted.bam; do
    SAMPLE_NAME=$(basename "$SAMPLE" .sorted.bam)

    # ヒトのリード長ごとのカウント取得
    samtools view "$SAMPLE" | awk '{print length($10)}' | sort | uniq -c | awk '{print $2, $1}' > "bam/${SAMPLE_NAME}_read_length.txt"

    # 各リード長ごとの正規化係数を計算
    awk 'NR==FNR {median[$1]=$2; next} {if ($1 in median) print $1, median[$1] / $2}' dm6_median_read_length.txt "bam/${SAMPLE_NAME}_read_length.txt" > "bam/${SAMPLE_NAME}_norm_factors.txt"

    # bedGraph の作成
    echo "Generating normalized bedGraph for $SAMPLE_NAME"
    bedtools genomecov -bg -ibam "$SAMPLE" -strand + -scale "$(awk '{sum += $2} END {print sum/NR}' bam/${SAMPLE_NAME}_norm_factors.txt)" -3 > "bw_readLength_normalized/${SAMPLE_NAME}_plus_norm.bedGraph"
    bedtools genomecov -bg -ibam "$SAMPLE" -strand - -scale "$(awk '{sum += $2} END {print sum/NR}' bam/${SAMPLE_NAME}_norm_factors.txt)" -3 > "bw_readLength_normalized/${SAMPLE_NAME}_minus_norm.bedGraph"

    # bigWig 変換
    bedGraphToBigWig "bw_readLength_normalized/${SAMPLE_NAME}_plus_norm.bedGraph" $CHROM_SIZES "bw_readLength_normalized/${SAMPLE_NAME}_plus_norm.bw"
    bedGraphToBigWig "bw_readLength_normalized/${SAMPLE_NAME}_minus_norm.bedGraph" $CHROM_SIZES "bw_readLength_normalized/${SAMPLE_NAME}_minus_norm.bw"

done

echo "Normalization completed!"