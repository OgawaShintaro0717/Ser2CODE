
cd ~/NGS/20240923si7SKpol2ChIP_GSE147055
cat SRR_Acc_List.txt | while read line;do cmd="fasterq-dump --split-files ${line}; gzip ${line}*fastq"; eval ${cmd}; done

cd ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastq
fasterq-dump --split-files SRR11313877
SRR11313858
SRR11313859
SRR11313860
SRR11313861
SRR11313862
SRR11313864
SRR11313866
SRR11313868
SRR11313870
SRR11313871
SRR11313876
SRR11313877

mkdir ~/NGS/20240923si7SKpol2ChIP_GSE147055 
mkdir fastqc
fastqc -o fastqc ~/NGS/20240923si7SKpol2ChIP_GSE147055/SRR11313858_1.fastq.gz ~/NGS/20240923si7SKpol2ChIP_GSE147055/SRR11313858_2.fastq.gz
fastqc -o fastqc fastq/*.fastq.gz
mkdir fastp
fastp --in1 fastq/SRR11313858_1.fastq.gz --in2 fastq/SRR11313858_2.fastq.gz --out1 fastp/SRR11313858_1.trim.fastq.gz --out2 fastp/SRR11313858_2.trim.fastq.gz --html fastp/SRR11313858.html 
fastp --in1 fastq/SRR11313859_1.fastq.gz --in2 fastq/SRR11313859_2.fastq.gz --out1 fastp/SRR11313859_1.trim.fastq.gz --out2 fastp/SRR11313859_2.trim.fastq.gz --html fastp/SRR11313859.html 
fastp --in1 fastq/SRR11313860_1.fastq.gz --in2 fastq/SRR11313860_2.fastq.gz --out1 fastp/SRR11313860_1.trim.fastq.gz --out2 fastp/SRR11313860_2.trim.fastq.gz --html fastp/SRR11313860.html 
fastp --in1 fastq/SRR11313861_1.fastq.gz --in2 fastq/SRR11313861_2.fastq.gz --out1 fastp/SRR11313861_1.trim.fastq.gz --out2 fastp/SRR11313861_2.trim.fastq.gz --html fastp/SRR11313861.html 
fastp --in1 fastq/SRR11313862_1.fastq.gz --in2 fastq/SRR11313862_2.fastq.gz --out1 fastp/SRR11313862_1.trim.fastq.gz --out2 fastp/SRR11313862_2.trim.fastq.gz --html fastp/SRR11313862.html 
fastp --in1 fastq/SRR11313864_1.fastq.gz --in2 fastq/SRR11313864_2.fastq.gz --out1 fastp/SRR11313864_1.trim.fastq.gz --out2 fastp/SRR11313864_2.trim.fastq.gz --html fastp/SRR11313864.html 
fastp --in1 fastq/SRR11313866_1.fastq.gz --in2 fastq/SRR11313866_2.fastq.gz --out1 fastp/SRR11313866_1.trim.fastq.gz --out2 fastp/SRR11313866_2.trim.fastq.gz --html fastp/SRR11313866.html 
fastp --in1 fastq/SRR11313868_1.fastq.gz --in2 fastq/SRR11313868_2.fastq.gz --out1 fastp/SRR11313868_1.trim.fastq.gz --out2 fastp/SRR11313868_2.trim.fastq.gz --html fastp/SRR11313868.html 
fastp --in1 fastq/SRR11313870_1.fastq.gz --in2 fastq/SRR11313870_2.fastq.gz --out1 fastp/SRR11313870_1.trim.fastq.gz --out2 fastp/SRR11313870_2.trim.fastq.gz --html fastp/SRR11313870.html 
fastp --in1 fastq/SRR11313871_1.fastq.gz --in2 fastq/SRR11313871_2.fastq.gz --out1 fastp/SRR11313871_1.trim.fastq.gz --out2 fastp/SRR11313871_2.trim.fastq.gz --html fastp/SRR11313871.html 
fastp --in1 fastq/SRR11313876_1.fastq.gz --in2 fastq/SRR11313876_2.fastq.gz --out1 fastp/SRR11313876_1.trim.fastq.gz --out2 fastp/SRR11313876_2.trim.fastq.gz --html fastp/SRR11313876.html 
fastp --in1 fastq/SRR11313877_1.fastq.gz --in2 fastq/SRR11313877_2.fastq.gz --out1 fastp/SRR11313877_1.trim.fastq.gz --out2 fastp/SRR11313877_2.trim.fastq.gz --html fastp/SRR11313877.html 

mkdir bowtie2
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313858_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313858_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313858.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313859_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313859_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313859.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313860_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313860_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313860.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313861_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313861_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313861.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313862_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313862_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313862.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313864_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313864_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313864.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313866_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313866_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313866.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313868_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313868_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313868.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313870_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313870_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313870.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313871_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313871_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313871.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313876_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313876_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313876.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313877_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313877_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim.sam

samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313858.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313858.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313858.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313859.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313859.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313859.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313860.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313860.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313860.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313861.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313861.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313861.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313862.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313862.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313862.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313864.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313864.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313864.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313866.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313866.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313866.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313868.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313868.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313868.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313870.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313870.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313870.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313871.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313871.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313871.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313876.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313876.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313876.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim.bam

samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313858.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313859.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313860.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313861.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313862.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313864.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313866.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313868.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313870.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313871.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313876.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim.bam

mkdir bigwig
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313858.trim.bam -o ./bigwig/SRR11313858.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313859.trim.bam -o ./bigwig/SRR11313859.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313860.trim.bam -o ./bigwig/SRR11313860.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313861.trim.bam -o ./bigwig/SRR11313861.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313862.trim.bam -o ./bigwig/SRR11313862.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313864.trim.bam -o ./bigwig/SRR11313864.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313866.trim.bam -o ./bigwig/SRR11313866.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313868.trim.bam -o ./bigwig/SRR11313868.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313870.trim.bam -o ./bigwig/SRR11313870.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313871.trim.bam -o ./bigwig/SRR11313871.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313876.trim.bam -o ./bigwig/SRR11313876.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313877.trim.bam -o ./bigwig/SRR11313877.trim.bw -of bigwig --normalizeUsing CPM



computeMatrix scale-regions -p 16 -S SRR5818147_7SKrep3Aligned.sortedByCoord.out.bw SRR5818143_CTrep3Aligned.sortedByCoord.out.bw -R hg38_majorRNU.bed --outFileName ~/NGS/20230106_Rn7SKseq/7SKRNAseq.txt.gz -a 1000 -b 1000 --skipZeros
plotProfile --yMin 0 -m ~/NGS/20230106_Rn7SKseq/7SKRNAseq.txt.gz --perGroup --plotTitle RNA-seq --samplesLabel si7SK siCT -out ~/NGS/20230106_Rn7SKseq/7SKRNAseq.txt.pdf


fastp --in1 fastq/SRR11313859_1.fastq.gz --in2 fastq/SRR11313859_2.fastq.gz --out1 fastp/SRR11313859_1.trim.fastq.gz --out2 fastp/SRR11313859_2.trim.fastq.gz --html fastp/SRR11313859.html 
fastp --in1 fastq/SRR11313860_1.fastq.gz --in2 fastq/SRR11313860_2.fastq.gz --out1 fastp/SRR11313860_1.trim.fastq.gz --out2 fastp/SRR11313860_2.trim.fastq.gz --html fastp/SRR11313860.html 
fastp --in1 fastq/SRR11313861_1.fastq.gz --in2 fastq/SRR11313861_2.fastq.gz --out1 fastp/SRR11313861_1.trim.fastq.gz --out2 fastp/SRR11313861_2.trim.fastq.gz --html fastp/SRR11313861.html 
fastp --in1 fastq/SRR11313862_1.fastq.gz --in2 fastq/SRR11313862_2.fastq.gz --out1 fastp/SRR11313862_1.trim.fastq.gz --out2 fastp/SRR11313862_2.trim.fastq.gz --html fastp/SRR11313862.html 
fastp --in1 fastq/SRR11313864_1.fastq.gz --in2 fastq/SRR11313864_2.fastq.gz --out1 fastp/SRR11313864_1.trim.fastq.gz --out2 fastp/SRR11313864_2.trim.fastq.gz --html fastp/SRR11313864.html 
fastp --in1 fastq/SRR11313866_1.fastq.gz --in2 fastq/SRR11313866_2.fastq.gz --out1 fastp/SRR11313866_1.trim.fastq.gz --out2 fastp/SRR11313866_2.trim.fastq.gz --html fastp/SRR11313866.html 

bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313859_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313859_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313859.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313860_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313860_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313860.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313861_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313861_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313861.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313862_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313862_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313862.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313864_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313864_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313864.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313866_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313866_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313866.trim.sam

samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313859.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313859.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313859.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313860.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313860.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313860.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313861.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313861.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313861.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313862.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313862.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313862.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313864.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313864.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313864.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313866.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313866.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313866.trim.bam

samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313859.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313860.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313861.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313862.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313864.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313866.trim.bam

bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313859.trim.bam -o ./bigwig/SRR11313859.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313860.trim.bam -o ./bigwig/SRR11313860.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313861.trim.bam -o ./bigwig/SRR11313861.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313862.trim.bam -o ./bigwig/SRR11313862.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313864.trim.bam -o ./bigwig/SRR11313864.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313866.trim.bam -o ./bigwig/SRR11313866.trim.bw -of bigwig --normalizeUsing CPM

fastp --in1 fastq/SRR11313868_1.fastq.gz --in2 fastq/SRR11313868_2.fastq.gz --out1 fastp/SRR11313868_1.trim.fastq.gz --out2 fastp/SRR11313868_2.trim.fastq.gz --html fastp/SRR11313868.html 
fastp --in1 fastq/SRR11313870_1.fastq.gz --in2 fastq/SRR11313870_2.fastq.gz --out1 fastp/SRR11313870_1.trim.fastq.gz --out2 fastp/SRR11313870_2.trim.fastq.gz --html fastp/SRR11313870.html 
fastp --in1 fastq/SRR11313871_1.fastq.gz --in2 fastq/SRR11313871_2.fastq.gz --out1 fastp/SRR11313871_1.trim.fastq.gz --out2 fastp/SRR11313871_2.trim.fastq.gz --html fastp/SRR11313871.html 
fastp --in1 fastq/SRR11313876_1.fastq.gz --in2 fastq/SRR11313876_2.fastq.gz --out1 fastp/SRR11313876_1.trim.fastq.gz --out2 fastp/SRR11313876_2.trim.fastq.gz --html fastp/SRR11313876.html 
fastp --in1 fastq/SRR11313877_1.fastq.gz --in2 fastq/SRR11313877_2.fastq.gz --out1 fastp/SRR11313877_1.trim.fastq.gz --out2 fastp/SRR11313877_2.trim.fastq.gz --html fastp/SRR11313877.html 

bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313868_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313868_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313868.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313870_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313870_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313870.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313871_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313871_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313871.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313876_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313876_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313876.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313877_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313877_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim.sam

samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313868.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313868.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313868.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313870.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313870.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313870.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313871.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313871.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313871.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313876.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313876.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313876.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim.bam

samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313868.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313870.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313871.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313876.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim.bam

bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313868.trim.bam -o ./bigwig/SRR11313868.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313870.trim.bam -o ./bigwig/SRR11313870.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313871.trim.bam -o ./bigwig/SRR11313871.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313876.trim.bam -o ./bigwig/SRR11313876.trim.bw -of bigwig --normalizeUsing CPM
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313877.trim.bam -o ./bigwig/SRR11313877.trim.bw -of bigwig --normalizeUsing CPM

fastp --in1 fastq/SRR11313877_1.fastq.gz --in2 fastq/SRR11313877_2.fastq.gz --out1 fastp/SRR11313877_1.trim.fastq.gz --out2 fastp/SRR11313877_2.trim.fastq.gz --html fastp/SRR11313877.html 
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -1 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313877_1.trim.fastq.gz -2 ~/NGS/20240923si7SKpol2ChIP_GSE147055/fastp/SRR11313877_2.trim.fastq.gz -S ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim.sam
samtools view -bhS -F 0x4 ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim.sam | samtools sort -T ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim - > ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim.bam
samtools index ~/NGS/20240923si7SKpol2ChIP_GSE147055/bowtie2/SRR11313877.trim.bam
bamCoverage -p 16 -bs 1 -b ./bowtie2/SRR11313877.trim.bam -o ./bigwig/SRR11313877.trim.bw -of bigwig --normalizeUsing CPM

cd ./bigwig
computeMatrix scale-regions -p 16 -S SRR11313870.trim.bw SRR11313871.trim.bw SRR11313876.trim.bw SRR11313877.trim.bw -R hg38_majorRNU.bed --outFileName pol2ChIP_seq.txt.gz -a 1000 -b 1000 --skipZeros
plotProfile --yMin 0 -m pol2ChIP_seq.txt.gz --perGroup --plotTitle "ChIPseq" --samplesLabel 7SK_1 7SK_2 CT_1 CT_2 -out pol2ChIP_seq.txt.pdf
plotProfile --yMin 0 --yMax 30 -m pol2ChIP_seq.txt.gz --perGroup --plotTitle "ChIPseq" --samplesLabel 7SK_1 7SK_2 CT_1 CT_2 -out pol2ChIP_seq_max30.txt.pdf

computeMatrix scale-regions -p 16 -S SRR11313871.trim.bw SRR11313876.trim.bw -R hg38_majorRNU.bed --outFileName pol2ChIP_seq_2.txt.gz -a 1000 -b 1000 --skipZeros
plotProfile --yMin 0 -m pol2ChIP_seq_2.txt.gz --perGroup --plotTitle "ChIPseq" --samplesLabel 7SK CT -out pol2ChIP_seq_2.txt.pdf
plotProfile --yMin 0 --yMax 30 -m pol2ChIP_seq_2.txt.gz --perGroup --plotTitle "ChIPseq" --samplesLabel 7SK CT -out pol2ChIP_seq_2_max30.txt.pdf

bigwigCompare -p 24 -b1 SRR11313876.trim.bw -b2 SRR11313877.trim.bw --skipNAs --operation mean --skipZeroOverZero -bs 1  -o CT_KO_average.bw
bigwigCompare -p 24 -b1 SRR11313870.trim.bw -b2 SRR11313871.trim.bw --skipNAs --operation mean --skipZeroOverZero -bs 1  -o 7SK_KO_average.bw

computeMatrix scale-regions -p 16 -S SRR11313860.trim.bw SRR11313861.trim.bw SRR11313866.trim.bw SRR11313868.trim.bw -R hg38_majorRNU.bed --outFileName 7SKChIRP_seq.txt.gz -a 1000 -b 1000 --skipZeros
plotProfile --yMin 0 -m 7SKChIRP_seq.txt.gz --perGroup --plotTitle "ChIRPseq" --samplesLabel input1 input2 7SK_1 7SK_2 -out 7SKChIRP_seq.txt.pdf

computeMatrix scale-regions -p 16 -S SRR11313860.trim.bw SRR11313866.trim.bw -R hg38_majorRNU.bed --outFileName 7SKChIRP_seq_1.txt.gz -a 1000 -b 1000 --skipZeros
plotProfile --yMin 0 -m 7SKChIRP_seq_1.txt.gz --perGroup --plotTitle "ChIRPseq" --samplesLabel input1 7SK_1 -out 7SKChIRP_seq_1.txt.pdf

computeMatrix scale-regions -p 16 -S SRR11313861.trim.bw SRR11313868.trim.bw -R hg38_majorRNU.bed --outFileName 7SKChIRP_seq_2.txt.gz -a 1000 -b 1000 --skipZeros
plotProfile --yMin 0 -m 7SKChIRP_seq_2.txt.gz --perGroup --plotTitle "ChIRPseq" --samplesLabel input2 7SK_2 -out 7SKChIRP_seq_2.txt.pdf


