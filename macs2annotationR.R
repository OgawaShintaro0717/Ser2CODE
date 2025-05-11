"C:\Documents\ICE2_ICE1_Pol2_Coilin\macs2\DMSO_ChIP_peaks.csv"
data <- read.csv("C:/Documents/ICE2_ICE1_Pol2_Coilin/macs2/DMSO_ChIP_peaks.csv")
head(data)
list.files("C:/Documents/ICE2_ICE1_Pol2_Coilin/macs2")
"C:\Documents\ICE2_ICE1_Pol2_Coilin\macs2\DMSO_ChIP_peaks.narrowPeak"
setwd("C:/Documents/ICE2_ICE1_Pol2_Coilin/macs2")
library(biomaRt)
library(ChIPseeker)
library(GenomicRanges)
library(clusterProfiler)
BiocManager::install("clusterProfiler")
a
package.version("ggplot2")

# BioconductorからTxDbオブジェクトをインストール
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
a
# TxDbパッケージをロード
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# narrowPeakファイルを読み込み
peak_file <- "C:/Documents/ICE2_ICE1_Pol2_Coilin/macs2/DMSO_ChIP_peaks.narrowPeak"
peaks <- read.table(peak_file, header = FALSE)
head(peaks)

# peakデータフレームをGenomicRangesオブジェクトに変換
peak_gr <- GRanges(seqnames = peaks$V1, 
                   ranges = IRanges(start = peaks$V2, end = peaks$V3),
                   peak_id = peaks$V4)
# ChIPseekerで遺伝子アノテーションを取得
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene  # hg38の場合（適宜変更）
peak_annotation <- annotatePeak(peak_gr, 
                                TxDb = txdb, 
                                tssRegion = c(-500, 500), 
                                annoDb = "org.Hs.eg.db")
# peak_annotationをデータフレームに変換
peak_annotation_df <- as.data.frame(peak_annotation)
# データフレームの先頭を確認
head(peak_annotation_df)
# 保存先のファイルパスを指定
output_file <- "C:/Documents/ICE2_ICE1_Pol2_Coilin/macs2/DMSO_ChIP_peaks.annotation.csv"
# CSVファイルとして保存
write.csv(peak_annotation_df, file = output_file, row.names = FALSE)

# 遺伝子タイプ（gene_biotype）を取得するためにbiomaRtを使用
# Ensembl BioMartに接続
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# 取得したい属性を指定 (gene_biotype)
attributes <- c("ensembl_gene_id", "gene_biotype")

# Ensemblの遺伝子IDに基づいて遺伝子タイプを取得
gene_biotype_info <- getBM(attributes = attributes, mart = ensembl)

# peak_annotation_df の geneId と gene_biotype_info の ensembl_gene_id をキーにして結合
peak_annotation_df <- merge(peak_annotation_df, gene_biotype_info, by.x = "ENSEMBL", by.y = "ensembl_gene_id", all.x = TRUE)

# 結果の確認
head(peak_annotation_df)

# 保存先のファイルパスを指定
output_file <- "C:/Documents/ICE2_ICE1_Pol2_Coilin/macs2/DMSO_ChIP_peaks.annotation_with_biotype.csv"

# CSVファイルとして保存
write.csv(peak_annotation_df, file = output_file, row.names = FALSE)

# dplyrパッケージをロード
library(dplyr)

# distanceToTSSが1000以下のものだけを抽出
filtered_peak_annotation_df <- peak_annotation_df %>%
  filter(abs(distanceToTSS) <= 1000)

# 結果の確認
head(filtered_peak_annotation_df)
########################################################################################################
#ﾊﾟｲﾁｬｰﾄ
library(ggplot2)
# gene_biotypeの頻度をカウント
gene_biotype_counts <- filtered_peak_annotation_df %>%
  count(gene_biotype) %>%
  arrange(desc(n))

# 各セクションの割合（%）を計算
gene_biotype_counts <- gene_biotype_counts %>%
  mutate(percentage = n / sum(n) * 100)
# 色の指定（青－緑系の色）
biotype_colors <- c(
  "protein_coding" = "#56B1F7",   # 青
  "lncRNA" = "#32A852",           # 緑
  "snoRNA" = "#1E90FF",           # 青系
  "snRNA" = "#00FA9A",            # 緑系
  "miRNA" = "#3CB371"             # ミドリ系
  # 必要に応じて他のタイプの色を追加してください
)

# パイチャートの作成
ggplot(gene_biotype_counts, aes(x = "", y = n, fill = gene_biotype)) +
  geom_bar(stat = "identity", width = 1) +  # バーの幅を1に設定して円形に
  coord_polar(theta = "y") +  # 極座標で円形に
  labs(title = "Gene Biotype Distribution") +
  theme_void() +  # 軸や背景を非表示
  theme(legend.title = element_blank()) +  # 凡例タイトルを非表示
  scale_fill_manual(values = biotype_colors) +  # 手動で色を設定
  geom_text(aes(label = paste0(round(percentage, 1), "%")),  # パーセント表示
            position = position_stack(vjust = 0.5),  # セクションの中央に配置
            color = "white", size = 5)
######################################################################################################
# `GENENAME` が重複している行を除去
# `gene_biotype` が NA の行も除去
cleaned_peak_annotation_df <- peak_annotation_df %>%
  filter(!duplicated(GENENAME)) %>%  # GENENAMEが重複している行を除去
  filter(!is.na(gene_biotype))  # gene_biotypeがNAの行を除去

# 結果を確認
head(cleaned_peak_annotation_df)

# gene_biotypeの頻度をカウント
gene_biotype_counts <- cleaned_peak_annotation_df %>%
  count(gene_biotype) %>%
  arrange(desc(n))

# 各セクションの割合（%）を計算
gene_biotype_counts <- gene_biotype_counts %>%
  mutate(percentage = n / sum(n) * 100)

# 色の指定（青－緑系の色）
biotype_colors <- c(
  "protein_coding" = "#56B1F7",   # 青
  "lncRNA" = "#32A852",           # 緑
  "snoRNA" = "#1E90FF",           # 青系
  "snRNA" = "#00FA9A",            # 緑系
  "miRNA" = "#3CB371"             # ミドリ系
  # 必要に応じて他のタイプの色を追加してください
)

# パイチャートの作成
ggplot(gene_biotype_counts, aes(x = "", y = n, fill = gene_biotype)) +
  geom_bar(stat = "identity", width = 1) +  # バーの幅を1に設定して円形に
  coord_polar(theta = "y") +  # 極座標で円形に
  labs(title = "Gene Biotype Distribution") +
  theme_void() +  # 軸や背景を非表示
  theme(legend.title = element_blank()) +  # 凡例タイトルを非表示
  scale_fill_manual(values = biotype_colors) +  # 手動で色を設定
  geom_text(aes(label = paste0(round(percentage, 1), "%")),  # パーセント表示
            position = position_stack(vjust = 0.5),  # セクションの中央に配置
            color = "white", size = 5)
