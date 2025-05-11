#DL fastq
cd ~/NGS/20240729ChIRCseq
fasterq-dump SRR10315849
fasterq-dump SRR16612767
fasterq-dump SRR16612766
fasterq-dump SRR10315841
fasterq-dump SRR10315839
fasterq-dump SRR10315845

gzip SRR10315849.fastq
gzip SRR16612767.fastq
gzip SRR16612766.fastq
gzip SRR10315841.fastq
gzip SRR10315839.fastq
gzip SRR10315845.fastq

