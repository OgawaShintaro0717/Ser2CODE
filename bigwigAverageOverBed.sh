
cd ~/NGS/S2_S5_Pol2ChIP
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed
sudo mv bigWigAverageOverBed /usr/local/bin/
bigWigAverageOverBed ~/NGS/S2_S5_Pol2ChIP/bw/SRR11440180.trim.bw all_genes.bed SRR11440180_output.tab
awk 'BEGIN {print "ENSG,Size,CoveredBases,Sum,Mean"} {print $1 "," $2 "," $3 "," $4 "," $5}' SRR11440180_output.tab > SRR11440180_output.csv

#!/bin/bash

# .bw ファイルのリスト
bw_files=("SRR11440180.trim.bw" "SRR11440181.trim.bw" "SRR11440186.trim.bw" "SRR11440187.trim.bw" "SRR11440190.trim.bw" "SRR11440191.trim.bw" "SRR11440194.trim.bw" "SRR11440195.trim.bw")

# 各ファイルについて処理を実行
for bw_file in "${bw_files[@]}"; do
    # bigWigAverageOverBedを実行して出力を .tab ファイルに保存
    bigWigAverageOverBed ~/NGS/S2_S5_Pol2ChIP/bw/$bw_file all_genes.bed ${bw_file%.bw}_output.tab

    # awk で CSV に変換
    awk 'BEGIN {print "ENSG,Size,CoveredBases,Sum,Mean"} {print $1 "," $2 "," $3 "," $4 "," $5}' ${bw_file%.bw}_output.tab > ${bw_file%.bw}_output.csv
done




