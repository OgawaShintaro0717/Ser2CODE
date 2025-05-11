#computeMatrix
cd ~/NGS/20240526_4sU_si7SK_hist
computeMatrix scale-regions -p 12 -S ~/NGS/20240526_4sU_si7SK_hist/SRR5818152Aligned.sortedByCoord.out.bw ~/NGS/20240526_4sU_si7SK_hist/SRR5818155Aligned.sortedByCoord.out.bw -R ~/NGS/20240526_4sU_si7SK_hist/221120_hist.bed --outFileName 4sUsi7SK_hist_1.txt.gz -a 500 -b 500 --skipZeros
computeMatrix scale-regions -p 12 -S ~/NGS/20240526_4sU_si7SK_hist/SRR5818152Aligned.sortedByCoord.out.bw ~/NGS/20240526_4sU_si7SK_hist/SRR5818155Aligned.sortedByCoord.out.bw -R ~/NGS/20240526_4sU_si7SK_hist/ProteinCodingGene_hg38.bed --outFileName 4sUsi7SK_proteinCodinggene_1.txt.gz -a 500 -b 500 --skipZeros

\\wsl.localhost\Ubuntu\home\shintaro\NGS\20240526_4sU_si7SK_hist
20240526_4sU_si7SK_hist
SRR5818152Aligned.sortedByCoord.out.bw
~/NGS/20240526_4sU_si7SK_hist/221120_hist.bed

#plotProfile
plotProfile --perGroup --samplesLabel "siCT" "si7SK" -m ~/NGS/20240526_4sU_si7SK_hist/4sUsi7SK_hist_1.txt.gz -out ~/NGS/20240526_4sU_si7SK_hist/4sUsi7SK_hist_1.pdf
plotProfile --perGroup --samplesLabel "siCT" "si7SK" -m ~/NGS/20240526_4sU_si7SK_hist/4sUsi7SK_proteinCodinggene_1.txt.gz -out ~/NGS/20240526_4sU_si7SK_hist/4sUsi7SK_proteinCoding_1.pdf

plotProfile --perGroup --samplesLabel "siCT" "si7SK" --plotType heatmap --kmeans 3 -m ~/NGS/20240526_4sU_si7SK_hist/4sUsi7SK_hist_1.txt.gz -out ~/NGS/20240526_4sU_si7SK_hist/4sUsi7SK_hist_heat_1.pdf
plotProfile --perGroup --samplesLabel "siCT" "si7SK" --plotType heatmap --kmeans 3 -m ~/NGS/20240526_4sU_si7SK_hist/4sUsi7SK_proteinCodinggene_1.txt.gz -out ~/NGS/20240526_4sU_si7SK_hist/4sUsi7SK_proteinCodinggene_heat_1.pdf



