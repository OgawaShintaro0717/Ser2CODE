#RNAseq
cd ~/NGS/20240626ICE1AFF4_RNAseq
/AFF4_rep1Aligned.toTranscriptome.out.bam"

~/STAR-2.7.11b/bin/Linux_x86_64/STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference --genomeFastaFiles ~/NGS/20240626ICE1AFF4_RNAseq/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile ~/NGS/20240626ICE1AFF4_RNAseq/Homo_sapiens.GRCh38.113.gtf
gunzip ~/NGS/20240626ICE1AFF4_RNAseq/fastq/AFF4_rep1_RNA_1.fq.gz ~/NGS/20240626ICE1AFF4_RNAseq/fastq/AFF4_rep1_RNA_2.fq.gz
cd ~/NGS/20240626ICE1AFF4_RNAseq/fastq
gunzip AFF4_rep2_RNA_1.fq.gz AFF4_rep2_RNA_2.fq.gz AFF4_rep3_RNA_1.fq.gz AFF4_rep3_RNA_2.fq.gz ICE1_rep1_RNA_1.fq.gz ICE1_rep1_RNA_2.fq.gz ICE1_rep2_RNA_1.fq.gz ICE1_rep2_RNA_2.fq.gz ICE1_rep3_RNA_1.fq.gz ICE1_rep3_RNA_2.fq.gz

#mapping
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/20240626ICE1AFF4_RNAseq/fastq/AFF4_rep2_RNA_1.fq ~/NGS/20240626ICE1AFF4_RNAseq/fastq/AFF4_rep2_RNA_2.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix AFF4rep2_50 --quantMode TranscriptomeSAM
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/20240626ICE1AFF4_RNAseq/fastq/AFF4_rep3_RNA_1.fq ~/NGS/20240626ICE1AFF4_RNAseq/fastq/AFF4_rep3_RNA_2.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix AFF4rep3_50 --quantMode TranscriptomeSAM
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/20240626ICE1AFF4_RNAseq/fastq/ICE1_rep1_RNA_1.fq ~/NGS/20240626ICE1AFF4_RNAseq/fastq/ICE1_rep1_RNA_2.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ICE1rep1_50 --quantMode TranscriptomeSAM
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/20240626ICE1AFF4_RNAseq/fastq/ICE1_rep2_RNA_1.fq ~/NGS/20240626ICE1AFF4_RNAseq/fastq/ICE1_rep2_RNA_2.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ICE1rep2_50 --quantMode TranscriptomeSAM
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/20240626ICE1AFF4_RNAseq/fastq/ICE1_rep3_RNA_1.fq ~/NGS/20240626ICE1AFF4_RNAseq/fastq/ICE1_rep3_RNA_2.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ICE1rep3_50 --quantMode TranscriptomeSAM

~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/20240626ICE1AFF4_RNAseq/AFF4rep2_50Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference AFF4_rep2
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/20240626ICE1AFF4_RNAseq/AFF4rep3_50Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference AFF4_rep3
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/20240626ICE1AFF4_RNAseq/ICE1rep1_50Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference ICE1_rep1
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/20240626ICE1AFF4_RNAseq/ICE1rep2_50Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference ICE1_rep2
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/20240626ICE1AFF4_RNAseq/ICE1rep3_50Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference ICE1_rep3

~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/20240626ICE1AFF4_RNAseq/fastq/AFF4_rep1_RNA_1.fq ~/NGS/20240626ICE1AFF4_RNAseq/fastq/AFF4_rep1_RNA_2.fq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix AFF4rep1_50 --quantMode TranscriptomeSAM
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/20240626ICE1AFF4_RNAseq/AFF4rep1_50Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference AFF4_rep1

mkdir ~/NGS/202411137SK_KD_RNAseq/bam
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818141_1.fastq ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818141_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818141 --quantMode TranscriptomeSAM
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818142_1.fastq ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818142_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818142 --quantMode TranscriptomeSAM
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818143_1.fastq ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818143_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818143 --quantMode TranscriptomeSAM
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818144_1.fastq ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818144_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818144 --quantMode TranscriptomeSAM
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818145_1.fastq ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818145_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818145 --quantMode TranscriptomeSAM
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818146_1.fastq ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818146_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818146 --quantMode TranscriptomeSAM
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818147_1.fastq ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818147_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818147 --quantMode TranscriptomeSAM
~/STAR-2.7.11b/bin/Linux_x86_64_static/STAR --runThreadN 16 --outFilterMultimapNmax 50 --runMode alignReads --genomeDir ~/NGS/20240626ICE1AFF4_RNAseq/STAR_reference -c --readFilesIn ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818148_1.fastq ~/NGS/202411137SK_KD_RNAseq/fastq/SRR5818148_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818148 --quantMode TranscriptomeSAM

mkdir ~/NGS/202411137SK_KD_RNAseq/rsem
cd ~/NGS/202411137SK_KD_RNAseq/rsem
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818141Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference SRR5818141
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818142Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference SRR5818142
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818143Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference SRR5818143
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818144Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference SRR5818144
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818145Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference SRR5818145
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818146Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference SRR5818146
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818147Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference SRR5818147
~/RSEM-1.3.3/rsem-calculate-expression -p 16 --paired-end --bam ~/NGS/202411137SK_KD_RNAseq/bam/SRR5818148Aligned.toTranscriptome.out.bam  ~/NGS/20240626ICE1AFF4_RNAseq/RSEM_reference/RSEM_reference SRR5818148



~/RSEM-1.3.3/rsem-prepare-reference -p 16 --gtf Homo_sapiens.GRCh38.113.gtf Homo_sapiens.GRCh38.dna.primary_assembly.fa RSEM_reference/RSEM_reference
"\\wsl.localhost\Ubuntu\home\shintaro\NGS\ref\RSEM_reference\RSEM_reference.stat"


~/Documents/expression/tools/RSEM-1.3.3/rsem-calculate-expression -p 12 --bam ~/Documents/expression/20230717ICE2KO/bam/SRR16701690Aligned.toTranscriptome.out.bam ~/Documents/expression/ref/RSEM_reference/RSEM_reference SRR16701690
"\\wsl.localhost\Ubuntu\home\shintaro\RSEM-1.3.3\boost"

setwd("C:/Documents/241114si7SKRNAseq")/SRR5818141.csv
~/Documents/241114si7SKRNAseq/SRR5818141.csv


