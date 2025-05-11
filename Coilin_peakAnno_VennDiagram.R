
library(biomaRt)
library(ChIPseeker)
library(GenomicRanges)
library(clusterProfiler)
BiocManager::install("clusterProfiler")
a
library(GenomicRanges)


# narrowPeakファイルを読み込み
peak_file <- "C:/Documents/20241129CoilinGFP_ChIPmacs2_output_1/Coilin_ChIP_peaks.narrowPeak"
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
                                tssRegion = c(-3000, 3000), 
                                annoDb = "org.Hs.eg.db")
# peak_annotationをデータフレームに変換
peak_annotation_df <- as.data.frame(peak_annotation)
# データフレームの先頭を確認
head(peak_annotation_df)
# 保存先のファイルパスを指定
output_file <- "C:/Documents/20241129CoilinGFP_ChIPmacs2_output_1/Coilin_ChIP_peak_annotations.csv"
# CSVファイルとして保存
write.csv(peak_annotation_df, file = output_file, row.names = FALSE)

coilin <- read.csv("C:/Documents/20241129CoilinGFP_ChIPmacs2_output_1/Coilin_ChIP_peak_annotations.csv")
selevted_genes <- read.csv("C:/Documents/20241129CoilinGFP_ChIPmacs2_output_1/selected_genes_info.csv")

head(coilin)
head(selevted_genes)

# NAを除去して、重複を取り除いたENSG IDのリストを作成
coilin_ensg <- unique(coilin$ENSEMBL[!is.na(coilin$ENSEMBL)])  # coilinのENSG
selected_ensg <- unique(selevted_genes$ENSG[!is.na(selevted_genes$ENSG)])  # selected_genesのENSG

# Venn図の作成
venn.diagram(
  x = list(Coilin = coilin_ensg, Selected = selected_ensg),  # セットを指定
  category.names = c("Coilin", "Selected Genes"),  # カテゴリ名を設定
  filename = "C:/Documents/20241129CoilinGFP_ChIPmacs2_output_1/venn_diagram.svg",  # 保存するファイル名
  imagetype = "svg",  # 画像形式
  height = 5,  # 画像の高さ
  width = 5,  # 画像の幅
  fill = c(4, 7),  # Coilinを左（色4）、Selected Genesを右（色7）
  lty = 1,  # 線のスタイル
  scaled = TRUE,  # 面積を要素数に比例させる
  cex = 2,  # 数字ラベルのフォントサイズ
  cat.pos = c(180, 0),  # Coilinを左（180度）、Selected Genesを右（0度）
  cat.dist = c(0.05, 0.05),  # カテゴリラベルと円の距離
  cat.cex = c(1.2, 1.2),  # カテゴリラベルのフォントサイズ
  rotation.degree = 0,  # 円の回転を設定しない
  subcat.pos = c(0.5, -0.5),  # 各円の配置を明示的に指定
  subcat.dist = 0.05  # カテゴリラベルの距離を調整
)


# Venn図の作成
venn.diagram(
  x = list(Coilin = coilin_ensg, Selected = selected_ensg),  # セットを指定
  category.names = c("Coilin", "Selected Genes"),  # カテゴリ名を設定
  filename = "C:/Documents/20241129CoilinGFP_ChIPmacs2_output_1/venn_diagram.svg",  # 保存するファイル名
  imagetype = "svg",  # 画像形式
  height = 5,  # 画像の高さ
  width = 5,  # 画像の幅
  fill = c(4, 7),  # 円の色
  lty = 1, # 線のスタイル
  scaled = TRUE,#の面積を要素数に比例させない
  cex = 2,  # 数字ラベルのフォントサイズ
  cat.pos = c(330, 30),  # カテゴリラベルの位置
  cat.dist = c(0.05, 0.05),  # カテゴリラベルと円の距離
  cat.cex = c(1.2, 1.2)  # カテゴリラベルのフォントサイズ
)

