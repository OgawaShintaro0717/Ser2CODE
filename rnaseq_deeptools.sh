"E:\molecularbiology\NGSwithGTune\20240626ICE1AFF4_RNAseq"

cd ~/NGS/20240626ICE1AFF4_RNAseq
computeMatrix scale-regions -p 14 -S  ~/NGS/20240626ICE1AFF4_RNAseq/AFF4_rep2Aligned.sortedByCoord.out.bw ~/NGS/20240626ICE1AFF4_RNAseq/ICE1_rep2Aligned.sortedByCoord.out.bw -R ~/NGS/20240626ICE1AFF4_RNAseq/hg38_majorRNU.bed --outFileName ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_5.txt.gz  -a 1000 -b 200 --skipZeros
plotProfile --regionsLabel RNUgenes --plotType lines --perGroup -m ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_5.txt.gz --yMin 0 --samplesLabel AFF4 ICE1  -out ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_5.pdf
plotProfile --regionsLabel RNUgenes --perGroup -m ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_5.txt.gz --yMax 5 --yMin 0 --samplesLabel AFF4 ICE1  -out ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_5_1.pdf
plotProfile --regionsLabel RNUgenes --plotType heatmap --perGroup -m ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_5.txt.gz --yMin 0 --samplesLabel AFF4 ICE1  -out ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_5_2.pdf
plotProfile --regionsLabel RNUgenes --plotType heatmap --perGroup -m ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_5.txt.gz --yMax 5 --yMin 0 --samplesLabel AFF4 ICE1  -out ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_5_3.pdf

~/NGS/20240626ICE1AFF4_RNAseq/AFF4_rep1Aligned.sortedByCoord.out.bw
~/NGS/20240626ICE1AFF4_RNAseq/AFF4_rep2Aligned.sortedByCoord.out.bw
~/NGS/20240626ICE1AFF4_RNAseq/AFF4_rep3Aligned.sortedByCoord.out.bw

~/NGS/20240626ICE1AFF4_RNAseq/ICE1_rep1Aligned.sortedByCoord.out.bw
~/NGS/20240626ICE1AFF4_RNAseq/ICE1_rep2Aligned.sortedByCoord.out.bw
~/NGS/20240626ICE1AFF4_RNAseq/ICE1_rep3Aligned.sortedByCoord.out.bw

~/NGS/20240626ICE1AFF4_RNAseq/hg38_majorRNU.bed
~/NGS/20240626ICE1AFF4_RNAseq/ProteinCodingGene_hg38.bed

computeMatrix scale-regions -p 20 -S  ~/NGS/20240626ICE1AFF4_RNAseq/AFF4_rep2Aligned.sortedByCoord.out.bw ~/NGS/20240626ICE1AFF4_RNAseq/ICE1_rep2Aligned.sortedByCoord.out.bw -R ~/NGS/20240626ICE1AFF4_RNAseq/ProteinCodingGene_hg38.bed --outFileName ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_6.txt.gz  -a 1000 -b 200 --skipZeros
plotProfile --regionsLabel ProteinCodingGenes --yMin 0 --plotType lines --perGroup -m ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_6.txt.gz --samplesLabel AFF4 ICE1  -out ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_6.pdf
plotProfile --regionsLabel ProteinCodingGenes --yMin 0 --plotType heatmap --perGroup -m ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_6.txt.gz --samplesLabel AFF4 ICE1  -out ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_7heat.pdf

plotProfile --kmeans 3 --regionsLabel ProteinCodingGenes --yMin 0 --plotType lines --perGroup -m ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_6.txt.gz --samplesLabel AFF4 ICE1  -out ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_6_kemeans.pdf
plotProfile --kmeans 3 --regionsLabel ProteinCodingGenes --yMin 0 --plotType heatmap --perGroup -m ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_6.txt.gz --samplesLabel AFF4 ICE1  -out ~/NGS/20240626ICE1AFF4_RNAseq/test_RNU_7heat_kemeans.pdf


cd ~/NGS/20250319CRC_RNU
computeMatrix scale-regions -p 16 -S  SRR13307408.bw SRR13307410.bw -R hg38_majorRNU.bed --outFileName CRC_RNU.txt.gz  -a 1000 -b 1000 --skipZeros
plotProfile --regionsLabel RNU_Genes --yMin 0 --plotType lines --perGroup -m CRC_RNU.txt.gz --samplesLabel NCM460 HCT116  -out NCM_CRC.pdf
plotProfile --regionsLabel RNU_Genes --yMin 0 --yMax 10 --plotType lines --perGroup -m CRC_RNU.txt.gz --samplesLabel NCM460 HCT116  -out NCM_CRC_max10.pdf




bigWigAverageOverBed SRR13307410.bw hg38_majorRNU_fixed.bed gene_body_counts_HCT.tab
bigWigAverageOverBed SRR13307410.bw hg38_majorRNU_TES_downstream_fixed.bed downstream_counts_HCT.tab

paste gene_body_counts_HCT.tab downstream_counts_HCT.tab | \
awk '{if ($5 > 0) print $1, $5, $10, $10/$5}' > read_through_ratio_HCT.tab
awk 'BEGIN {OFS=","} {print $1, $2, $3, $4}' read_through_ratio.tab > read_through_ratio.csv
awk 'BEGIN {OFS=","} {print $1, $2, $3, $4}' read_through_ratio_HCT.tab > read_through_ratio_HCT.csv
