"\\wsl.localhost\Ubuntu\home\shintaro\NGS\20211217_7SK_ChIRP\202410057SK_ChIRP_GSE69141_SRP058600.sh"
cd ~/NGS/20211217_7SK_ChIRP
mkdir fastq
cd ~/NGS/20211217_7SK_ChIRP/fastq
cat ~/NGS/20211217_7SK_ChIRP/SRR_Acc_List.txt | while read line;do cmd="fasterq-dump ${line}; gzip ${line}*fastq"; eval ${cmd}; done
fasterq-dump SRR2036989
gzip SRR2036989.fastq

mkdir ~/NGS/20211217_7SK_ChIRP 
mkdir fastqc
fastqc -o fastqc fastq/*.fastq.gz
SRR2036981
SRR2036982
SRR2036984
SRR2036983
SRR2036985
SRR2036986
SRR2036987
SRR2036988
SRR2036989

81-84 mES
83-86 hES
87-89 Hela
even odd input


mkdir ~/NGS/20211217_7SK_ChIRP/fastp
cd ~/NGS/20211217_7SK_ChIRP/fastp
cat ~/NGS/20211217_7SK_ChIRP/SRR_Acc_List.txt | while read line;do cmd="fastp -i ${line}.fastq.gz -o ~/NGS/20211217_7SK_ChIRP/fastp/${line}.trim.fastq.gz --html ~/NGS/20211217_7SK_ChIRP/fastp/${line}.html "; eval ${cmd}; done
fastp -i fastq/SRR2036989.fastq.gz -o fastp/SRR2036989.trim.fastq.gz --html fastp/SRR2036989.html 

mkdir ~/NGS/20211217_7SK_ChIRP/bowtie2
cat ~/NGS/20211217_7SK_ChIRP/SRR_Acc_List_mm.txt | while read line;do cmd="bowtie2 -p 16 -x ~/bowtie2_index/GRCm39/GRCm39 -U ~/NGS/20211217_7SK_ChIRP/fastp/${line}.trim.fastq.gz -S ~/NGS/20211217_7SK_ChIRP/bowtie2/${line}.trim.sam "; eval ${cmd}; done
cat ~/NGS/20211217_7SK_ChIRP/SRR_Acc_List_hs.txt | while read line;do cmd="bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -U ~/NGS/20211217_7SK_ChIRP/fastp/${line}.trim.fastq.gz -S ~/NGS/20211217_7SK_ChIRP/bowtie2/${line}.trim.sam "; eval ${cmd}; done
bowtie2 -p 16 -x ~/bowtie2_index/GRCm39/GRCm39 -U ~/NGS/20211217_7SK_ChIRP/fastp/SRR2036981.fastq.trim.gz -S ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR11313858.trim.sam
bowtie2 -p 16 -x ~/bowtie2_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -U ~/NGS/20211217_7SK_ChIRP/fastp/SRR2036989.fastq.trim.gz -S ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036989.trim.sam


mkdir ~/NGS/20211217_7SK_ChIRP/bowtie2
cd ~/NGS/20211217_7SK_ChIRP/bowtie2
cat ~/NGS/20211217_7SK_ChIRP/SRR_Acc_List_hs.txt | while read line;do cmd="samtools view -bhS -F 0x4 ~/NGS/20211217_7SK_ChIRP/bowtie2/${line}.trim.sam ; samtools sort -T ~/NGS/20211217_7SK_ChIRP/bowtie2/${line}.trim - > ~/NGS/bowtie2/${line}.trim.bam "; eval ${cmd}; done

samtools view -bhS -F 0x4 ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036981.trim.sam | samtools sort -T ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036981.trim - > ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036981.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036982.trim.sam | samtools sort -T ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036982.trim - > ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036982.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036983.trim.sam | samtools sort -T ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036983.trim - > ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036983.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036984.trim.sam | samtools sort -T ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036984.trim - > ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036984.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036985.trim.sam | samtools sort -T ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036985.trim - > ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036985.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036986.trim.sam | samtools sort -T ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036986.trim - > ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036986.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036987.trim.sam | samtools sort -T ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036987.trim - > ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036987.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036988.trim.sam | samtools sort -T ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036988.trim - > ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036988.trim.bam
samtools view -bhS -F 0x4 ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036989.trim.sam | samtools sort -T ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036989.trim - > ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036989.trim.bam

cat ~/NGS/20211217_7SK_ChIRP/SRR_Acc_List.txt | while read line;do cmd="samtools index ~/NGS/20211217_7SK_ChIRP/bowtie2/${line}.trim.bam "; eval ${cmd}; done
samtools index ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036989.trim.bam

mkdir ~/NGS/20211217_7SK_ChIRP/bigwig
cat ~/NGS/20211217_7SK_ChIRP/SRR_Acc_List_hs.txt | while read line;do cmd="bamCoverage -p 16 -bs 1 -b ~/NGS/20211217_7SK_ChIRP/bowtie2/${line}.trim.bam -o ~/NGS/20211217_7SK_ChIRP/bigwig/${line}.trim.bw -of bigwig --normalizeUsing CPM "; eval ${cmd}; done

mkdir ~/NGS/20211217_7SK_ChIRP/macs2
source ~/tools/MACS2/bin/activate 
macs2 callpeak -t ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036987.trim.bam -c ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036989.trim.bam -f BAM -g hs -n 7SK_ChIRP_HeLa_even --outdir macs2 -B -q 0.01
macs2 callpeak -t ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036988.trim.bam -c ~/NGS/20211217_7SK_ChIRP/bowtie2/SRR2036989.trim.bam -f BAM -g hs -n 7SK_ChIRP_HeLa_odd --outdir macs2 -B -q 0.01


computeMatrix scale-regions -p 16 -S SRR5818147_7SKrep3Aligned.sortedByCoord.out.bw SRR5818143_CTrep3Aligned.sortedByCoord.out.bw -R hg38_majorRNU.bed --outFileName ~/NGS/20230106_Rn7SKseq/7SKRNAseq.txt.gz -a 1000 -b 1000 --skipZeros
plotProfile --yMin 0 -m ~/NGS/20230106_Rn7SKseq/7SKRNAseq.txt.gz --perGroup --plotTitle RNA-seq --samplesLabel si7SK siCT -out ~/NGS/20230106_Rn7SKseq/7SKRNAseq.txt.pdf

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


