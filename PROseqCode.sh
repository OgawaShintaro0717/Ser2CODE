
########################################################
#!/bin/bash

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
            -S output/${SAMPLE}.sam --threads $THREADS --very-sensitive --no-mixed --no-discordant

    ## 3. BAM処理（ソート & 重複除去）
    samtools view -bS output/${SAMPLE}.sam | samtools sort -o bam/${SAMPLE}.sorted.bam
    samtools index bam/${SAMPLE}.sorted.bam
    rm output/${SAMPLE}.sam  # SAM削除してディスク節約

    ## 4. 5'末端のカバレッジ計算（bedGraph → bigWig）
    bedtools genomecov -bg -ibam bam/${SAMPLE}.sorted.bam -strand + -3 > bw/${SAMPLE}_3plus.bedGraph
    bedtools genomecov -bg -ibam bam/${SAMPLE}.sorted.bam -strand - -3 > bw/${SAMPLE}_3minus.bedGraph

    # bigWig変換
    bedGraphToBigWig bw/${SAMPLE}_3plus.bedGraph $CHROM_SIZES bw/${SAMPLE}_3plus.bw
    bedGraphToBigWig bw/${SAMPLE}_3minus.bedGraph $CHROM_SIZES bw/${SAMPLE}_3minus.bw

done

echo "Pipeline completed!"

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

    ## 4. 5'末端のカバレッジ計算（bedGraph → bigWig）
    bedtools genomecov -bg -ibam bam/${SAMPLE}.sorted.bam -strand + -3 > bw/${SAMPLE}_3plus.bedGraph
    bedtools genomecov -bg -ibam bam/${SAMPLE}.sorted.bam -strand - -3 > bw/${SAMPLE}_3minus.bedGraph

    # bigWig変換
    bedGraphToBigWig bw/${SAMPLE}_3plus.bedGraph $CHROM_SIZES bw/${SAMPLE}_3plus.bw
    bedGraphToBigWig bw/${SAMPLE}_3minus.bedGraph $CHROM_SIZES bw/${SAMPLE}_3minus.bw

done

echo "Pipeline completed!"



bedGraphToBigWig ~/NGS/250219PROseq/bw/1_ICE1_1_S1_L001_plus.bedGraph ~/NGS/250219PROseq/hg38.chrom.sizes ~/NGS/250219PROseq/bw/1_ICE1_1_S1_L001_plus.bw
bedGraphToBigWig ~/NGS/250219PROseq/bw/1_ICE1_1_S1_L001_minus.bedGraph ~/NGS/250219PROseq/hg38.chrom.sizes ~/NGS/250219PROseq/bw/1_ICE1_1_S1_L001_minus.bw

awk '{ if ($5 == "-") $4 = -$4; print }' /home/shintaro/NGS/250219PROseq/bw/1_ICE1_1_S1_L001_minus.bedGraph > /home/shintaro/NGS/250219PROseq/bw/1_ICE1_1_S1_L001_minus_reflected.bedGraph
bedGraphToBigWig /home/shintaro/NGS/250219PROseq/bw/1_ICE1_1_S1_L001_minus_reflected.bedGraph /home/shintaro/NGS/250219PROseq/hg38.chrom.sizes /home/shintaro/NGS/250219PROseq/bw/1_ICE1_1_S1_L001_minus_reflected.bw

cd ~/NGS/250219PROseq/bw_readLength_normalized
awk '!($4 ~ /^RNU6/ || $4 ~ /^RNU7/)' hg38_snRNA_filtered_plus_only.bed > hg38_snRNA_filtered_plus_RNU67filtered_bed_file.bed


# マイナスストランド用のメタジーン
computeMatrix scale-regions -S /path/to/minus.bw -R /path/to/hg38_majorRNU.bed --upstream 0 --downstream 1000 -o /path/to/minus_matrix.gz

# マイナスストランドをプラスに揃える
computeMatrix scale-regions -S /path/to/minus.bw -R /path/to/hg38_majorRNU.bed --upstream 1000 --downstream 0 --flip -o /path/to/minus_flipped_matrix.gz

# プラスストランドとマイナスストランドを合わせてメタジーンプロットを描く
plotProfile --perGroup -m plus_matrix.gz -o plus_profile.pdf --yMin 0

plotProfile -m /path/to/minus_flipped_matrix.gz -t "Minus Strand Reads" -o /path/to/minus_profile.png

plotProfile -m /path/to/plus_matrix.gz /path/to/minus_flipped_matrix.gz \
            -t "Plus and Minus Strand Reads" \
            --samplesLabel "Plus Strand" "Minus Strand" \
            -o /path/to/combined_profile.png

plotHeatmap -m /path/to/plus_matrix.gz /path/to/minus_flipped_matrix.gz \
            --colorMap RdBu \
            --samplesLabel "Plus Strand" "Minus Strand" \
            -o /path/to/combined_heatmap.png

computeMatrix reference-point \
  -S /path/to/plus.bw /path/to/minus.bw \
  -R /path/to/hg38_majorRNU.bed \
  --beforeRegionStartLength 1000 \
  --afterRegionStartLength 1000 \
  --referencePoint TSS \
  -o /path/to/output_matrix.gz

cd ~/bowtie2_index/BDGP6
ls ~/bowtie2_index/BDGP6

#!/bin/bash

# 変数定義
BOWTIE_INDEX=~/bowtie2_index/BDGP6/BDGP6  # ハエのBowtie2インデックス
CHROM_SIZES=dm6.chrom.sizes  # ハエのchrom.sizes
THREADS=16  # 並列実行スレッド数


# サンプルごとに処理
for R1 in ~/NGS/250219PROseq/*_R1_001.fastq.gz; do
    SAMPLE=$(basename $R1 _R1_001.fastq.gz)
    R2=~/NGS/250219PROseq/${SAMPLE}_R2_001.fastq.gz

    echo "Processing: $SAMPLE (mapping to dm6)"

    ## 2. アライメント (bowtie2)
    bowtie2 -x $BOWTIE_INDEX \
            -1 trimmed_logs_dm6/${SAMPLE}_R1_trimmed.fastq.gz \
            -2 trimmed_logs_dm6/${SAMPLE}_R2_trimmed.fastq.gz \
            -S output_dm6/${SAMPLE}.sam \
            --threads $THREADS --very-sensitive --no-mixed --no-discordant

    ## 3. BAM処理（ソート & 重複除去）
    samtools view -bS output_dm6/${SAMPLE}.sam | samtools sort -o bam_dm6/${SAMPLE}.sorted.bam
    samtools index bam_dm6/${SAMPLE}.sorted.bam
    rm output_dm6/${SAMPLE}.sam  # SAM削除してディスク節約

    ## 4. 5'末端のカバレッジ計算（bedGraph → bigWig）
    bedtools genomecov -bg -ibam bam_dm6/${SAMPLE}.sorted.bam -strand + -5 > bw_dm6/${SAMPLE}_plus.bedGraph
    bedtools genomecov -bg -ibam bam_dm6/${SAMPLE}.sorted.bam -strand - -5 > bw_dm6/${SAMPLE}_minus.bedGraph

    # bigWig変換
    bedGraphToBigWig bw_dm6/${SAMPLE}_plus.bedGraph $CHROM_SIZES bw_dm6/${SAMPLE}_plus.bw
    bedGraphToBigWig bw_dm6/${SAMPLE}_minus.bedGraph $CHROM_SIZES bw_dm6/${SAMPLE}_minus.bw

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
#上記のコードを参考に、bam_noRepeatsの中にあるヒトのbamファイルに対し、bam_dm6_noRepeatsの中にあるハエのbamファイルを用いて正規化をかけたい。出力としてはbwとbedGraphでいきたい
#一応コードの意味としてはリードの長さごとの数を反映して正規化をかけている感じ




cd ~/NGS/250219PROseq/bam_dm6
mkdir -p bam_dm_noRepeats

for BAM in *.sorted.bam; do
    SAMPLE=$(basename "$BAM" .sorted.bam)
    FILTERED_BAM="bam_dm_noRepeats/${SAMPLE}_noRepeats.bam"
    
    echo "Filtering $SAMPLE..."
    samtools view -b -q 30 -F 0x100 -o "$FILTERED_BAM" "$BAM"

    echo "Indexing $FILTERED_BAM..."
    samtools index "$FILTERED_BAM"
done


cd ~/NGS/bedfile
# GTFファイルからprotein_codingの遺伝子を抽出し、BED形式に変換
awk '$3 == "gene" && $0 ~ /gene_type "protein_coding"/ {split($9, a, ";"); for(i in a) if(a[i] ~ /gene_name/) gene_name=a[i]; print $1 "\t" $4-1 "\t" $5 "\t" gene_name "\t" 0 "\t" $7}' ~/NGS/bedfile/gencode.v45.annotation.gtf > protein_coding_genes.bed
awk '$5 == "+" {print $0}' protein_coding_genes.bed > protein_coding_genes_plus.bed
awk '$5 == "-" {print $0}' protein_coding_genes.bed > protein_coding_genes_minus.bed

head ~/NGS/bedfile/gencode.v45.annotation.gtf

cd ~/NGS/250219PROseq/bw_readLength_normalized
# プラスストランド用のメタジーン
computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R hg38_majorRNUplus.bed --upstream 1000 --downstream 1000 --skipZeros -o plus_matrix.gz
plotProfile --perGroup -m plus_matrix.gz -o RNU_Max500_plus_profile.pdf --yMin 0 --yMax 500
plotProfile --perGroup -m plus_matrix.gz -o RNU_plus_profile.pdf --yMin 0

computeMatrix scale-regions -S 1_ICE1_1_S1_L001_plus_norm.bw 4_ICE1AFF4_2_S4_L001_plus_norm.bw -R protein_coding_genes_top500.bed --upstream 1000 --downstream 1000 --skipZeros -o Protein_plus_matrix.gz
#plotProfile --perGroup -m Protein_plus_matrix.gz -o Protein_Max500_plus_profile.pdf --yMin 0 --yMax 500
plotProfile --perGroup -m Protein_plus_matrix.gz -o Protein_plus_profile.pdf --yMin 0 


ls ~/NGS/250219PROseq/bam
head ~/NGS/bedfile/gencode.v45.annotation.gtf
head ~/NGS/bedfile/hg38_majorRNUplus.bed
gtf_file = "/home/your_username/NGS/bedfile/gencode.v45.annotation.gtf"
realpath /home/your_username/NGS/bedfile/gencode.v45.annotation.gtf
ls -l /home/your_username/NGS/bedfile/gencode.v45.annotation.gtf

import pandas as pd
import pysam
import re

# GTFファイルのパス
gtf_file = "/home/shintaro/NGS/bedfile/gencode.v45.annotation.gtf"


# BAMファイルリスト
bam_files = [
    "1_ICE1_1_S1_L001.sorted.bam",
    "2_ICE1_2_S2_L001.sorted.bam",
    "3_ICE1AFF4_1_S3_L001.sorted.bam",
    "4_ICE1AFF4_2_S4_L001.sorted.bam"
]

# BAMファイルのディレクトリ
bam_dir = "/home/shintaro/250219PROseq/bam"

# BEDフォーマットの出力ファイル
output_bed = "protein_coding_genes_sorted.bed"

# 1. Protein Coding Gene の抽出
protein_coding_genes = []
with open(gtf_file, "r") as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if fields[2] == "gene":
            match = re.search(r'gene_type "(.*?)"', fields[8])
            if match and match.group(1) == "protein_coding":
                gene_name = re.search(r'gene_name "(.*?)"', fields[8]).group(1)
                chrom, start, end, strand = fields[0], int(fields[3])-1, int(fields[4]), fields[6]
                protein_coding_genes.append([chrom, start, end, gene_name, strand])

# データフレーム化
genes_df = pd.DataFrame(protein_coding_genes, columns=["chr", "start", "end", "gene_name", "strand"])

# 2. BAMファイルからリード数カウント
gene_counts = {gene[3]: 0 for gene in protein_coding_genes}  # 初期化

for bam_file in bam_files:
    bam_path = f"{bam_dir}/{bam_file}"
    bam = pysam.AlignmentFile(bam_path, "rb")

    for gene in protein_coding_genes:
        chrom, start, end, gene_name, strand = gene
        count = bam.count(contig=chrom, start=start, end=end)
        gene_counts[gene_name] += count

    bam.close()

# 3. リード数順にソート
sorted_genes = sorted(protein_coding_genes, key=lambda x: gene_counts[x[3]], reverse=True)

# 4. BEDフォーマットで出力
with open(output_bed, "w") as bed:
    for gene in sorted_genes:
        chrom, start, end, gene_name, strand = gene
        bed.write(f"{chrom}\t{start}\t{end}\t0\t{gene_name}\t{strand}\n")

print(f"BEDファイルが作成されました: {output_bed}")


computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R protein_coding_plus_strand_top100.bed --upstream 1000 --downstream 1000 --skipZeros -o Protein_plus_matrix.gz
plotProfile --perGroup -m Protein_plus_matrix.gz -o Protein_plus_profile.pdf --yMin 0 --samplesLabel "AFF4" "ICE1"

computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R snRNA_plus_strand_top100.bed --upstream 1000 --downstream 1000 --skipZeros -o RNU_plus_matrix.gz
plotProfile --perGroup -m RNU_plus_matrix.gz -o RNU_plus_profile.pdf --yMin 0 --samplesLabel "AFF4" "ICE1"

computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R protein_coding_genes_plus.bed --upstream 1000 --downstream 1000 --skipZeros -o Protein_plus_all_matrix.gz
plotProfile --perGroup -m Protein_plus_all_matrix.gz -o Protein_plus_all_profile.pdf --yMin 0 --samplesLabel "AFF4" "ICE1"

computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R final_clean_genes_top1000.bed --upstream 1000 --downstream 1000 --skipZeros -o final_clean_genes_top1000_matrix.gz
plotProfile --perGroup -m final_clean_genes_top1000_matrix.gz -o final_clean_genes_top1000.pdf --yMin 0 --samplesLabel "AFF4" "ICE1"

head protein_coding_plus_strand_top1000.bed 

awk '{
    if ($6 == "+") {
        start = ($2 - 1000 < 0) ? 0 : $2 - 1000;
        end = $2;
    } else {
        start = $3;
        end = $3 + 1000;
    }
    print $1"\t"start"\t"end"\t"$5"\t.\t"$6;
}'  protein_coding_genes_plus_with_names.bed > upstream_1kb_with_names.bed
head upstream_1kb.bed
awk '{
    if ($6 == "+") {
        start = $3;
        end = $3 + 1000;
    } else {
        start = ($2 - 1000 < 0) ? 0 : $2 - 1000;
        end = $2;
    }
    print $1"\t"start"\t"end"\t"$5"\t.\t"$6;
}'  protein_coding_genes_plus_with_names.bed > downstream_1kb_with_names.bed

awk '$3 == "gene" {
    match($0, /gene_name "([^"]+)"/, arr);
    print $1"\t"$4-1"\t"$5"\t"arr[1]"\t.\t"$7;
}' ~/NGS/bedfile/gencode.v45.annotation.gtf > all_genes.bed

bedtools intersect -v -a upstream_1kb_with_names.bed -b all_genes.bed > upstream_clean.bed
head upstream_clean.bed
bedtools intersect -v -a downstream_1kb_with_names.bed -b all_genes.bed > downstream_clean.bed

cut -f4 upstream_clean.bed | sort > upstream_genes.txt
cut -f4 downstream_clean.bed | sort > downstream_genes.txt

comm -12 upstream_genes.txt downstream_genes.txt > genes_clean_updown.txt

grep -Ff genes_clean_updown.txt protein_coding_genes_plus.bed > final_clean_genes_updown.bed
head final_clean_genes_updown.bed

computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R final_clean_genes_top1000.bed --upstream 1000 --downstream 1000 --skipZeros -o final_clean_genes_top1000_matrix.gz
plotProfile --perGroup -m final_clean_genes_top1000_matrix.gz -o final_clean_genes_top1000_1.pdf --yMin 1.9 --yMax 5.5 --samplesLabel "AFF4" "ICE1"


(base) shintaro@Double-OO:~/NGS/250219PROseq/bw_readLength_normalized$ head protein_coding_genes_plus_with_names.bed
chr1    65419   71585   OR4F5   .       +
chr1    923923  944575  SAMD11  .       +
chr1    960584  965719  KLHL17  .       +
chr1    966482  975865  PLEKHN1 .       +
chr1    1001138 1014540 ISG15   .       +
chr1    1020120 1056118 AGRN    .       +
chr1    1173880 1197936 TTLL10  .       +
chr1    1232237 1235041 B3GALT6 .       +
chr1    1280436 1292029 SCNN1D  .       +
chr1    1308597 1311677 PUSL1   .       +
(base) shintaro@Double-OO:~/NGS/250219PROseq/bw_readLength_normalized$ head upstream_1kb.bed
chr1    64419   65419   OR4F5   .       +
chr1    922923  923923  SAMD11  .       +
chr1    959584  960584  KLHL17  .       +
chr1    965482  966482  PLEKHN1 .       +
chr1    1000138 1001138 ISG15   .       +
chr1    1019120 1020120 AGRN    .       +
chr1    1172880 1173880 TTLL10  .       +
chr1    1231237 1232237 B3GALT6 .       +
chr1    1279436 1280436 SCNN1D  .       +
chr1    1307597 1308597 PUSL1   .       +
(base) shintaro@Double-OO:~/NGS/250219PROseq/bw_readLength_normalized$ head downstream_1kb.bed
chr1    71585   72585   OR4F5   .       +
chr1    944575  945575  SAMD11  .       +
chr1    965719  966719  KLHL17  .       +
chr1    975865  976865  PLEKHN1 .       +
chr1    1014540 1015540 ISG15   .       +
chr1    1056118 1057118 AGRN    .       +
chr1    1197936 1198936 TTLL10  .       +
chr1    1235041 1236041 B3GALT6 .       +
chr1    1292029 1293029 SCNN1D  .       +
chr1    1311677 1312677 PUSL1   .       +
(base) shintaro@Double-OO:~/NGS/250219PROseq/bw_readLength_normalized$ head all_genes.bed
chr1    11868   14409   DDX11L2 .       +
chr1    12009   13670   DDX11L1 .       +
chr1    14695   24886   WASH7P  .       -
chr1    17368   17436   MIR6859-1       .       -
chr1    29553   31109   MIR1302-2HG     .       +
chr1    30365   30503   MIR1302-2       .       +
chr1    34553   36081   FAM138A .       -
chr1    52472   53312   OR4G4P  .       +
chr1    57597   64116   ENSG00000290826 .       +
chr1    62948   63887   OR4G11P .       +
####################################################################################################
awk '{
    if ($6 == "+") {
        start = ($2 - 1000 < 0) ? 0 : $2 - 1000;
        end = $2;
    } else {
        start = $3;
        end = $3 + 1000;
    }
    print $1"\t"start"\t"end"\t"$5"\t.\t"$6;
}' protein_coding_plus_strand_top1000_noChrM.bed > upstream_1kb_regions_top1000.bed

awk '{
    if ($6 == "+") {
        start = $3;
        end = $3 + 1000;
    } else {
        start = ($2 - 1000 < 0) ? 0 : $2 - 1000;
        end = $2;
    }
    print $1"\t"start"\t"end"\t"$5"\t.\t"$6;
}'  protein_coding_plus_strand_top1000_noChrM.bed > downstream_1kb_regions_top1000.bed

awk '$3 == "gene" {
    match($0, /gene_name "([^"]+)"/, arr);
    print $1"\t"$4-1"\t"$5"\t"arr[1]"\t.\t"$7;
}' ~/NGS/bedfile/gencode.v45.annotation.gtf > all_genes.bed

bedtools intersect -v -a upstream_1kb_regions_top1000.bed -b all_genes.bed > upstream_clean_top1000.bed
bedtools intersect -v -a downstream_1kb_regions_top1000.bed -b all_genes.bed > downstream_clean_top1000.bed


cut -f4 upstream_clean_top1000.bed > gene_names_no_upstream_overlap_top1000.txt
cut -f4 downstream_clean_top1000.bed > gene_names_no_downstream_overlap_top1000.txt
comm -12 gene_names_no_upstream_overlap_top1000.txt gene_names_no_downstream_overlap_top1000.txt > top1000_genes_clean_updown.txt

grep -Ff top1000_genes_clean_updown.txt protein_coding_plus_strand_top1000_noChrM.bed > final_clean_genes_top1000.bed

computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R final_clean_genes_top1000.bed --upstream 1000 --downstream 1000 --skipZeros -o final_clean_genes_top1000_matrix_2.gz
plotProfile --perGroup -m final_clean_genes_top1000_matrix_2.gz -o final_clean_genes_top1000_2.pdf --yMin 0 --samplesLabel "AFF4" "ICE1"
###################################################################################################
awk '$1 != "chrM"' protein_coding_plus_strand_top1000.bed > protein_coding_plus_strand_top1000_noChrM.bed
awk '{
    if ($6 == "+") {
        start = ($2 - 1000 < 0) ? 0 : $2 - 1000;
        end = $2;
    } else {
        start = $3;
        end = $3 + 1000;
    }
    print $1"\t"start"\t"end"\t"$5"\t.\t"$6;
}' protein_coding_plus_strand_top1000_noChrM.bed > upstream_1kb_regions_top1000.bed
awk '$3 == "gene" {
    match($0, /gene_name "([^"]+)"/, arr);
    print $1"\t"$4-1"\t"$5"\t"arr[1]"\t.\t"$7;
}' ~/NGS/bedfile/gencode.v45.annotation.gtf > all_genes.bed
bedtools intersect -v -a upstream_1kb_regions_top1000.bed -b all_genes.bed > upstream_clean_top1000.bed

cut -f4 upstream_clean_top1000.bed > gene_names_no_upstream_overlap_top1000.txt

grep -Ff gene_names_no_upstream_overlap_top1000.txt protein_coding_plus_strand_top1000_noChrM.bed > final_clean_genes_top1000.bed

computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R final_clean_genes_top1000.bed --upstream 1000 --downstream 1000 --skipZeros -o final_clean_genes_top1000_matrix_3.gz
plotProfile --perGroup -m final_clean_genes_top1000_matrix_3.gz -o final_clean_genes_top1000_3.pdf --yMin 1.9 --yMax 5.5 --samplesLabel "AFF4" "ICE1"
####################################################################################################

awk '$1 != "chrM"' protein_coding_plus_strand_top500.bed > protein_coding_plus_strand_top500_noChrM.bed
awk '{
    if ($6 == "+") {
        start = ($2 - 1000 < 0) ? 0 : $2 - 1000;
        end = $2;
    } else {
        start = $3;
        end = $3 + 1000;
    }
    print $1"\t"start"\t"end"\t"$5"\t.\t"$6;
}' protein_coding_plus_strand_top500_noChrM.bed > upstream_1kb_regions_top500.bed

bedtools intersect -v -a upstream_1kb_regions_top500.bed -b all_genes.bed > upstream_clean_top500.bed

cut -f4 upstream_clean_top500.bed > gene_names_no_upstream_overlap_top500.txt

grep -Ff gene_names_no_upstream_overlap_top500.txt protein_coding_plus_strand_top500_noChrM.bed > final_clean_genes_top500.bed

computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R final_clean_genes_top500.bed --upstream 1000 --downstream 1000 --skipZeros -o final_clean_genes_top500_matrix_3.gz
plotProfile --perGroup -m final_clean_genes_top500_matrix_3.gz -o final_clean_genes_top500_3.pdf --yMin 1.9  --samplesLabel "AFF4" "ICE1"

computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R final_clean_genes_top500.bed --upstream 500 --downstream 500 --skipZeros -o final_clean_genes_top500_matrix_4.gz
plotProfile --perGroup -m final_clean_genes_top500_matrix_4.gz -o final_clean_genes_top500_4.pdf --yMin 1.9  --samplesLabel "AFF4" "ICE1"

####################################################################################################

awk '$1 != "chrM"' protein_coding_plus_strand_top500.bed > protein_coding_plus_strand_top500_noChrM.bed
awk '{
    if ($6 == "+") {
        start = ($2 - 2000 < 0) ? 0 : $2 - 1000;
        end = $2;
    } else {
        start = $3;
        end = $3 + 2000;
    }
    print $1"\t"start"\t"end"\t"$5"\t.\t"$6;
}' protein_coding_plus_strand_top500_noChrM.bed > upstream_2kb_regions_top500.bed

bedtools intersect -v -a upstream_2kb_regions_top500.bed -b all_genes.bed > upstream_2kb_clean_top500.bed

cut -f4 upstream_2kb_clean_top500.bed > gene_names_2kb_no_upstream_overlap_top500.txt

grep -Ff gene_names_2kb_no_upstream_overlap_top500.txt protein_coding_plus_strand_top500_noChrM.bed > final_clean_genes_2kb_top500.bed

computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R final_clean_genes_2kb_top500.bed --upstream 1000 --downstream 1000 --skipZeros -o final_clean_genes_2kb_top500_matrix_3.gz
plotProfile --perGroup -m final_clean_genes_2kb_top500_matrix_3.gz -o final_clean_genes_2kb_top500_3.pdf --yMin 1.9  --samplesLabel "AFF4" "ICE1"

computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R final_clean_genes_top500.bed --upstream 500 --downstream 500 --skipZeros -o final_clean_genes_top500_matrix_4.gz
plotProfile --perGroup -m final_clean_genes_top500_matrix_4.gz -o final_clean_genes_top500_4.pdf --yMin 1.9  --samplesLabel "AFF4" "ICE1"

####################################################################################################
head protein_coding_genes_sorted.bed
awk '$6 == "+"' protein_coding_genes_sorted.bed > protein_coding_genes_sorted_plus.bed
head protein_coding_genes_sorted_plus.bed

awk '$1 != "chrM"'  protein_coding_genes_sorted_plus.bed >  protein_coding_genes_sorted_plus.bed_noChrM.bed

awk '{
    if ($6 == "+") {
        start = ($2 - 1000 < 0) ? 0 : $2 - 1000;
        end = $2;
    } else {
        start = $3;
        end = $3 + 1000;
    }
    print $1"\t"start"\t"end"\t"$5"\t.\t"$6;
}' protein_coding_genes_sorted_plus.bed_noChrM.bed > upstream_1kb_regions_sorted.bed

awk '{
    if ($6 == "+") {
        start = $3;
        end = $3 + 1000;
    } else {
        start = ($2 - 1000 < 0) ? 0 : $2 - 1000;
        end = $2;
    }
    print $1"\t"start"\t"end"\t"$5"\t.\t"$6;
}'  protein_coding_genes_sorted_plus.bed_noChrM.bed > downstream_1kb_regions_sorted.bed


bedtools intersect -v -a upstream_1kb_regions_sorted.bed -b all_genes.bed > upstream_clean_sorted.bed
bedtools intersect -v -a downstream_1kb_regions_sorted.bed -b all_genes.bed > downstream_clean_sorted_.bed
head upstream_clean_sorted.bed
head downstream_clean_sorted_.bed
cut -f4 upstream_clean_sorted.bed > gene_names_no_upstream_overlap_sorted.txt
cut -f4 downstream_clean_sorted_.bed > gene_names_no_downstream_overlap_sorted.txt
head gene_names_no_upstream_overlap_sorted.txt
head gene_names_no_downstream_overlap_sorted.txt
comm -12 gene_names_no_upstream_overlap_sorted.txt gene_names_no_downstream_overlap_sorted.txt > sorted_genes_clean_updown.txt
head sorted_genes_clean_updown.txt
grep -Ff sorted_genes_clean_updown.txt protein_coding_genes_sorted_plus.bed_noChrM.bed > final_clean_genes_sorted.bed
head final_clean_genes_sorted.bed
computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R final_clean_genes_sorted.bed --upstream 1000 --downstream 1000 --skipZeros -o final_clean_genes_sorted_matrix_2.gz
plotProfile --perGroup -m final_clean_genes_sorted_matrix_2.gz -o final_clean_genes_sorted_2.pdf --yMin 0 --samplesLabel "AFF4" "ICE1"

##########################################################################################################################
awk '{
    if ($6 == "+") {
        start = ($2 - 1000 < 0) ? 0 : $2 - 1000;
        end = $2;
    } else {
        start = $3;
        end = $3 + 1000;
    }
    print $1"\t"start"\t"end"\t"$5"\t.\t"$6;
}' protein_coding_genes_sorted_plus.bed > upstream_1kb_regions_sorted_plus.bed

bedtools intersect -v -a upstream_1kb_regions_sorted_plus.bed -b all_genes.bed > upstream_clean_sorted_plus.bed

cut -f4 upstream_clean_sorted_plus.bed > gene_names_no_upstream_overlap_sorted_plus.txt
head gene_names_no_upstream_overlap_sorted_plus.txt
grep -Ff gene_names_no_upstream_overlap_sorted_plus.txt protein_coding_genes_sorted_plus.bed > final_clean_genes_sorted_plus.bed
head 
computeMatrix scale-regions -S 4_ICE1AFF4_2_S4_L001_plus_norm.bw 1_ICE1_1_S1_L001_plus_norm.bw -R final_clean_genes_sorted_plus.bed --upstream 1000 --downstream 1000 --skipZeros -o final_clean_genes_sorted_plus_matrix_3.gz
plotProfile --perGroup -m final_clean_genes_sorted_plus_matrix_3.gz -o final_clean_genes_sorted_plus_3.pdf --yMin 1.9  --samplesLabel "AFF4" "ICE1"


head 4_ICE1AFF4_2_S4_L001_minus_norm.bedGraph
awk '{print $1, $2, $3, -$4}' OFS='\t' 4_ICE1AFF4_2_S4_L001_minus_norm.bedGraph > 4_ICE1AFF4_2_S4_L001_minusInverted.bedGraph
bedGraphToBigWig 4_ICE1AFF4_2_S4_L001_minusInverted.bedGraph hg38.chrom.sizes 4_ICE1AFF4_2_S4_L001_minusInverted.bw

