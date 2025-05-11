if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "dplyr", "ggplot2", "Rtsne", "umap"))
a
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)
library(Rtsne)
library(umap)
library(DESeq2)
# データのクエリ作成
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# データのダウンロード
GDCdownload(query)
getwd()  # 現在の作業ディレクトリを確認
list.files("GDCdata")  # GDCdataフォルダがあるか確認
setwd("C:/Documents/20250201TCGAanalysis")
# データの準備
tcga_data <- GDCprepare(query)
tar -xvzf Sat_Feb__1_12_49_45_2025_0.tar.gz
tar -xvzf Sat_Feb__1_12_49_45_2025_1.tar.gz
tar -xvzf Sat_Feb__1_12_49_45_2025_2.tar.gz


data_A1 <- read.csv("C:/Documents/20241117RNAseqIntergration/ICE1AFF4rnaseqSTAR_/AFF4rep1_counts.csv")
head(data_A1)

tcga_counts <- assay(tcga_data)  # SummarizedExperimentオブジェクトからカウント行列を抽出
colnames(tcga_counts) <- gsub("\\..*", "", colnames(tcga_counts))  # サンプルIDを簡略化

# TCGAのサンプル情報を取得
tcga_sample_info <- colData(tcga_data)

# TCGAのサンプルをがん（Tumor）と正常（Normal）で分類
tcga_sample_info$condition <- ifelse(grepl("01A", tcga_sample_info$barcode), "Tumor", "Normal")
rownames(tcga_sample_info) <- colnames(tcga_counts)

# 必要なパッケージ
library(dplyr)

# データがあるディレクトリを指定
data_dir <- "C:/Documents/20241117RNAseqIntergration/ICE1AFF4rnaseqSTAR_"

# `.txt` ファイルのみ取得
ice1_aff4_files <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)

# カウントデータを読み込んで統合
count_list <- lapply(ice1_aff4_files, function(file) {
  df <- read.delim(file, row.names = 1, header = FALSE, sep = "\t")  # 遺伝子IDを行名に
  colnames(df) <- gsub(".txt", "", basename(file))  # ファイル名をサンプル名に
  return(df)
})

# すべてのデータをカラム方向に統合
ice1_aff4_counts <- do.call(cbind, count_list)

# 列名の確認（ファイル名が列名として適切に付いているか）
colnames(ice1_aff4_counts)

# ICE1/AFF4のサンプル情報を作成
ice1_aff4_sample_info <- data.frame(
  sample = colnames(ice1_aff4_counts),
  condition = ifelse(grepl("ICE1", colnames(ice1_aff4_counts)), "ICE1", "AFF4")  # ICE1とAFF4を判別
)
rownames(ice1_aff4_sample_info) <- colnames(ice1_aff4_counts)

# 統合データを保存
write.csv(ice1_aff4_counts, "C:/Documents/20241117RNAseqIntergration/ICE1_AFF4_counts_merged.csv")

# サンプル情報も保存
write.csv(ice1_aff4_sample_info, "C:/Documents/20241117RNAseqIntergration/ICE1_AFF4_sample_info.csv")

# データの確認
dim(ice1_aff4_counts)  # 遺伝子数 × サンプル数
head(ice1_aff4_counts)  # 最初の数行を表示


# ICE1/AFF4カウントデータを読み込み
ice1_aff4_counts <- read.csv("C:/Documents/20241117RNAseqIntergration/ICE1_AFF4_counts_merged.csv", row.names = 1)

# TCGAカウントデータを読み込み
tcga_counts <- assay(tcga_data)  # SummarizedExperimentオブジェクトからカウントデータを抽出
colnames(tcga_counts) <- gsub("\\..*", "", colnames(tcga_counts))  # TCGAのサンプル名を簡略化

# 共通遺伝子のみを取得
common_genes <- intersect(rownames(tcga_counts), rownames(ice1_aff4_counts))

# 共通遺伝子に絞る
tcga_counts <- tcga_counts[common_genes, ]
ice1_aff4_counts <- ice1_aff4_counts[common_genes, ]

# 統合データを作成
merged_counts <- cbind(tcga_counts, ice1_aff4_counts)

# サンプル情報を作成
tcga_sample_info <- data.frame(
  sample = colnames(tcga_counts),
  condition = ifelse(grepl("01A", colnames(tcga_counts)), "Tumor", "Normal")  # TCGAデータの条件（がん or 正常）
)

ice1_aff4_sample_info <- read.csv("C:/Documents/20241117RNAseqIntergration/ICE1_AFF4_sample_info.csv", row.names = 1)

merged_sample_info <- rbind(tcga_sample_info, ice1_aff4_sample_info)
rownames(merged_sample_info) <- merged_sample_info$sample
merged_sample_info <- merged_sample_info[,-1, drop = FALSE]  # 余計なカラムを削除

# 統合データを保存
write.csv(merged_counts, "C:/Documents/20241117RNAseqIntergration/Merged_TCGA_ICE1_AFF4_counts.csv")
write.csv(merged_sample_info, "C:/Documents/20241117RNAseqIntergration/Merged_TCGA_ICE1_AFF4_sample_info.csv")

# DESeq2データセットを作成
dds <- DESeqDataSetFromMatrix(countData = merged_counts, colData = merged_sample_info, design = ~ condition)

# 正規化
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# 正規化データを保存
write.csv(normalized_counts, "C:/Documents/20241117RNAseqIntergration/Normalized_TCGA_ICE1_AFF4_counts.csv")
###################################################################################################
library(ggplot2)

# PCA解析
vsd <- vst(dds, blind = TRUE)  # Variance Stabilizing Transformation（正規化 & 変動の安定化）
pca_data <- prcomp(t(assay(vsd)))

# データフレーム化
pca_df <- data.frame(
  PC1 = pca_data$x[,1],
  PC2 = pca_data$x[,2],
  condition = merged_sample_info$condition
)

# プロット
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA plot of TCGA and ICE1/AFF4 samples")
# UMAP解析
umap_res <- umap(t(assay(vsd)))

# UMAPデータフレーム化
umap_df <- data.frame(
  UMAP1 = umap_res$layout[,1],
  UMAP2 = umap_res$layout[,2],
  condition = merged_sample_info$condition
)

# プロット
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "UMAP plot of TCGA and ICE1/AFF4 samples")
######################################################################################################

# TCGAデータだけを抽出
tcga_only_counts <- tcga_counts  # ICE1/AFF4 を除外
tcga_only_sample_info <- tcga_sample_info  # ICE1/AFF4 のメタデータも除外

# DESeq2データセットを作成
dds_tcga <- DESeqDataSetFromMatrix(countData = tcga_only_counts, colData = tcga_only_sample_info, design = ~ condition)

# 正規化
dds_tcga <- DESeq(dds_tcga)
normalized_counts_tcga <- counts(dds_tcga, normalized = TRUE)

# 保存（もし再利用するなら）
write.csv(normalized_counts_tcga, "C:/Documents/20241117RNAseqIntergration/Normalized_TCGA_counts.csv")

# UMAP解析（正規化データを使用）
vsd_tcga <- vst(dds_tcga, blind = TRUE)  # Variance Stabilizing Transformation（正規化＆変動安定化）
umap_res_tcga <- umap(t(assay(vsd_tcga)))

# UMAPのデータフレームを作成
umap_df_tcga <- data.frame(
  UMAP1 = umap_res_tcga$layout[,1],
  UMAP2 = umap_res_tcga$layout[,2],
  condition = tcga_only_sample_info$condition
)

# UMAPプロット
ggplot(umap_df_tcga, aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "UMAP plot of TCGA-COAD (Tumor vs Normal)")
# UMAPの結果をCSVとして保存
write.csv(umap_df_tcga, "C:/Documents/20241117RNAseqIntergration/UMAP_TCGA_results.csv", row.names = FALSE)

# ICE1の発現量を取得（遺伝子名を確認して適宜変更）
ice1_gene <- "ENSG00000164151"  # 正確な遺伝子名に変更（Ensembl IDの場合は `ENSGXXXXXXX`）
ice1_expression <- normalized_counts_tcga[ice1_gene, , drop = FALSE]

# UMAPデータフレームにICE1の発現量を追加
umap_df_tcga$ICE1_expression <- as.numeric(ice1_expression[colnames(normalized_counts_tcga)])

#rownames(normalized_counts_tcga)[1:10]  # 先頭10個の遺伝子IDを確認
# バージョン番号（.XX）を削除
rownames(normalized_counts_tcga) <- gsub("\\..*", "", rownames(normalized_counts_tcga))


# ICE1の発現データを取得
if (ice1_gene %in% rownames(normalized_counts_tcga)) {
  ice1_expression <- normalized_counts_tcga[ice1_gene, , drop = FALSE]
} else {
  print("指定した遺伝子IDがデータ内に見つかりません")
}

# UMAPデータフレームに ICE1 の発現量を追加
umap_df_tcga$ICE1_expression <- as.numeric(ice1_expression[colnames(normalized_counts_tcga)])


umap_df_tcga$ICE1_expression <- as.numeric(ice1_expression[colnames(normalized_counts_tcga)])

# `NA` を 0 に置き換え（必要なら）
umap_df_tcga$ICE1_expression[is.na(umap_df_tcga$ICE1_expression)] <- 0


# ICE1 の発現データを取得（`rownames(umap_df_tcga)` の順番に並べる）
ice1_expression <- as.numeric(normalized_counts_tcga["ENSG00000164151", rownames(umap_df_tcga)])
# `NA` を 0 に置き換え（発現データがないサンプルのため）
ice1_expression[is.na(ice1_expression)] <- 0

# UMAPデータフレームに追加
umap_df_tcga$ICE1_expression <- ice1_expression

ggplot(umap_df_tcga, aes(x = UMAP1, y = UMAP2, color = ICE1_expression)) +
  geom_point(size = 3) +
  scale_color_gradientn(
    colors = c("red", "orange", "yellow", "lightblue", "blue"),  # 赤→青のグラデーション
    values = scales::rescale(c(min(umap_df_tcga$ICE1_expression), 1000, 1250, 1500, max(umap_df_tcga$ICE1_expression)))
  ) +
  theme_minimal() +
  labs(title = "UMAP plot colored by ICE1 expression",
       color = "ICE1 expression")
ggsave("UMAP_ICE1_expression.pdf", 
       plot = last_plot(),  # 直前に作成したプロットを保存
       width = 7, height = 5, dpi = 300, device = cairo_pdf)
###############################################################################################
# AFF4 の発現データを取得
aff4_gene <- "ENSG00000072364"

# UMAPのサンプル順でAFF4の発現データを取得
aff4_expression <- as.numeric(normalized_counts_tcga[aff4_gene, rownames(umap_df_tcga)])

# `NA` を 0 に置き換え（発現データがないサンプルのため）
aff4_expression[is.na(aff4_expression)] <- 0

# UMAPデータフレームに追加
umap_df_tcga$AFF4_expression <- aff4_expression
ggplot(umap_df_tcga, aes(x = UMAP1, y = UMAP2, color = AFF4_expression)) +
  geom_point(size = 3) +
  scale_color_gradientn(
    colors = c("red", "orange", "yellow", "lightblue", "blue"),  # 赤→青のグラデーション
    values = scales::rescale(c(min(umap_df_tcga$AFF4_expression), 1000, 2000, 3000, max(umap_df_tcga$AFF4_expression)))
  ) +
  theme_minimal() +
  labs(title = "UMAP plot colored by AFF4 expression",
       color = "AFF4 expression")
ggsave("UMAP_AFF4_expression.pdf", 
       plot = last_plot(),  # 直前に作成したプロットを保存
       width = 7, height = 5, dpi = 300, device = cairo_pdf)
#############################################################################################
# 発現比 r = AFF4 / ICE1 を計算（ゼロ除算を防ぐために +1 を加える）
umap_df_tcga$AFF4_ICE1_ratio <- (umap_df_tcga$AFF4_expression + 1) / (umap_df_tcga$ICE1_expression + 1)

# `log10` でスケールを調整（比が極端にならないように）
umap_df_tcga$AFF4_ICE1_ratio_log <- log10(umap_df_tcga$AFF4_ICE1_ratio)

# 発現比の統計情報を確認
summary(umap_df_tcga$AFF4_ICE1_ratio_log)
library(ggplot2)

ggplot(umap_df_tcga, aes(x = UMAP1, y = UMAP2, color = AFF4_ICE1_ratio_log)) +
  geom_point(size = 3) +
  scale_color_gradientn(
    colors = c("blue", "lightblue", "yellow", "orange", "red"),  # 小さい方（青）→ 大きい方（赤）
    values = scales::rescale(c(min(umap_df_tcga$AFF4_ICE1_ratio_log),  # 最小値（青）
                               0.2,  # 1st Qu.（ライトブルー）
                               0.3,  # 中央値（黄）
                               0.4,  # 3rd Qu.（オレンジ）
                               max(umap_df_tcga$AFF4_ICE1_ratio_log))) # 最大値（赤）
  ) +
  theme_minimal() +
  labs(title = "UMAP plot colored by AFF4 / ICE1 expression ratio (log10)",
       color = "log10(AFF4 / ICE1)")
ggsave("UMAP_AFF4_ICE1_ratio_red_high_blue_low.pdf", 
       plot = last_plot(),  
       width = 7, height = 5, dpi = 300, device = cairo_pdf)
##################################################################################################
#CDK9/ICE1
# CDK9 の発現データを取得
cdk9_gene <- "ENSG00000136807"

# UMAPのサンプル順でCDK9の発現データを取得
cdk9_expression <- as.numeric(normalized_counts_tcga[cdk9_gene, rownames(umap_df_tcga)])

# `NA` を 0 に置き換え（発現データがないサンプルのため）
cdk9_expression[is.na(cdk9_expression)] <- 0

# 発現比 r = CDK9 / ICE1 を計算（ゼロ除算を防ぐために +1 を加える）
umap_df_tcga$CDK9_ICE1_ratio <- (cdk9_expression + 1) / (umap_df_tcga$ICE1_expression + 1)

# `log10` でスケールを調整（比が極端にならないように）
umap_df_tcga$CDK9_ICE1_ratio_log <- log10(umap_df_tcga$CDK9_ICE1_ratio)

# 発現比の統計情報を確認
summary(umap_df_tcga$CDK9_ICE1_ratio_log)
library(ggplot2)

ggplot(umap_df_tcga, aes(x = UMAP1, y = UMAP2, color = CDK9_ICE1_ratio_log)) +
  geom_point(size = 3) +
  scale_color_gradientn(
    colors = c("blue", "lightblue", "yellow", "orange", "red"),  # 小さい方（青）→ 大きい方（赤）
    values = scales::rescale(c(min(umap_df_tcga$CDK9_ICE1_ratio_log),  
                               quantile(umap_df_tcga$CDK9_ICE1_ratio_log, 0.25),  # 1st Qu.（ライトブルー）
                               median(umap_df_tcga$CDK9_ICE1_ratio_log),  # 中央値（黄）
                               quantile(umap_df_tcga$CDK9_ICE1_ratio_log, 0.75),  # 3rd Qu.（オレンジ）
                               max(umap_df_tcga$CDK9_ICE1_ratio_log)))  # 最大値（赤）
  ) +
  theme_minimal() +
  labs(title = "UMAP plot colored by CDK9 / ICE1 expression ratio (log10)",
       color = "log10(CDK9 / ICE1)")
ggsave("UMAP_CDK9_ICE1_ratio_red_high_blue_low.pdf", 
       plot = last_plot(),  
       width = 7, height = 5, dpi = 300, device = cairo_pdf)
range(umap_df_tcga$CDK9_ICE1_ratio_log)
quantile(umap_df_tcga$CDK9_ICE1_ratio_log, probs = c(0.25, 0.5, 0.75))
mean(umap_df_tcga$CDK9_ICE1_ratio_log, na.rm = TRUE)
sd(umap_df_tcga$CDK9_ICE1_ratio_log, na.rm = TRUE)

#####################################################################################################
# 上位 10% を求めるための閾値を計算
threshold_top10 <- quantile(umap_df_tcga$CDK9_ICE1_ratio_log, 0.90)

# 上位 10% のサンプルを赤、それ以外を灰色にする
umap_df_tcga$CDK9_ICE1_highlight <- ifelse(umap_df_tcga$CDK9_ICE1_ratio_log >= threshold_top10, "red", "gray")

# 確認
table(umap_df_tcga$CDK9_ICE1_highlight)
ggplot(umap_df_tcga, aes(x = UMAP1, y = UMAP2, color = CDK9_ICE1_highlight)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red" = "red", "gray" = "gray")) +
  theme_minimal() +
  labs(title = "UMAP plot with top 10% CDK9 / ICE1 ratio highlighted",
       color = "Top 10% CDK9 / ICE1")

# 上位 10% の閾値を計算
threshold_top10 <- quantile(umap_df_tcga$CDK9_ICE1_ratio_log, 0.90)

# 上位 10% のサンプルを赤、それ以外を灰色にする
umap_df_tcga$CDK9_ICE1_highlight <- ifelse(umap_df_tcga$CDK9_ICE1_ratio_log >= threshold_top10, "red", "gray")

# 赤（上位10%）が後にプロットされるように並び替え
umap_df_tcga <- umap_df_tcga[order(umap_df_tcga$CDK9_ICE1_highlight == "red"), ]
ggplot(umap_df_tcga, aes(x = UMAP1, y = UMAP2, color = CDK9_ICE1_highlight)) +
  geom_point(size = 3) + 
  scale_color_manual(values = c("gray" = "gray", "red" = "red")) +  # 灰色→赤の順番を固定
  theme_minimal() +
  labs(title = "UMAP plot with top 10% CDK9 / ICE1 ratio highlighted",
       color = "CDK Group")
head(umap_df_tcga)
ggsave("UMAP_CDK9_ICE1_ratio_red_10high_gray_low.pdf", 
       plot = last_plot(),  
       width = 7, height = 5, dpi = 300, device = cairo_pdf)



# 癌の範囲をフィルタリング
cancer_region <- umap_df_tcga[
  umap_df_tcga$UMAP1 > -3 & umap_df_tcga$UMAP1 < 3 & 
    umap_df_tcga$UMAP2 > -5 & umap_df_tcga$UMAP2 < 2.5, ]

# 正常の範囲をフィルタリング
normal_region <- umap_df_tcga[
  umap_df_tcga$UMAP1 > 7.5 & umap_df_tcga$UMAP2 > 14, ]

# CDK9 / ICE1 の比率の上位 10% の閾値を計算
threshold_top10 <- quantile(umap_df_tcga$CDK9_ICE1_ratio_log, 0.90, na.rm = TRUE)

# 各グループ内で上位 10% の割合を計算
cancer_top10_ratio <- mean(cancer_region$CDK9_ICE1_ratio_log >= threshold_top10, na.rm = TRUE)
normal_top10_ratio <- mean(normal_region$CDK9_ICE1_ratio_log >= threshold_top10, na.rm = TRUE)

# 結果を表示
cancer_top10_ratio
normal_top10_ratio
# データフレームとして保存
top10_ratios <- data.frame(
  Group = c("Cancer", "Normal"),
  Top10_Ratio = c(cancer_top10_ratio, normal_top10_ratio)
)

# CSV に保存
write.csv(top10_ratios, "CDK9_ICE1_top10_ratios.csv", row.names = FALSE)

# 結果の確認
print(top10_ratios)
head(umap_df)
####################################################################################################
# CEA の Ensembl ID
cea_gene <- "ENSG00000170956"

# UMAPのサンプル順でceaの発現データを取得
cea_expression <- as.numeric(normalized_counts_tcga[cea_gene, rownames(umap_df_tcga)])

# `NA` を 0 に置き換え（発現データがないサンプルのため）
cea_expression[is.na(cea_expression)] <- 0

# UMAPデータフレームに CEA の発現データを追加
umap_df_tcga$CEA_expression <- cea_expression

# 発現値の統計情報を確認
summary(umap_df_tcga$CEA_expression)

ggplot(umap_df_tcga, aes(x = UMAP1, y = UMAP2, color = CEA_expression)) +
  geom_point(size = 3) +
  scale_color_gradientn(
    colors = c("blue", "lightblue", "yellow", "orange", "red"),  # 低発現（青）→ 高発現（赤）
    values = scales::rescale(c(min(umap_df_tcga$CEA_expression),  
                               quantile(umap_df_tcga$CEA_expression, 0.25),  # 1st Qu.（ライトブルー）
                               median(umap_df_tcga$CEA_expression),  # 中央値（黄）
                               quantile(umap_df_tcga$CEA_expression, 0.75),  # 3rd Qu.（オレンジ）
                               max(umap_df_tcga$CEA_expression)))  # 最大値（赤）
  ) +
  theme_minimal() +
  labs(title = "UMAP plot colored by CEA (ENSG00000170956) expression",
       color = "CEA Expression")
ggsave("UMAP_CEA_expression_70956.pdf", 
       plot = last_plot(),  
       width = 7, height = 5, dpi = 300, device = cairo_pdf)

ggplot(umap_df_tcga, aes(x = CEA_expression, y = CDK9_ICE1_ratio_log)) +
  geom_point(alpha = 0.6, color = "blue") +
  theme_minimal() +
  labs(title = "Scatter plot of CEA expression vs. CDK9 / ICE1 ratio",
       x = "CEA Expression",
       y = "log10(CDK9 / ICE1)") +
  geom_smooth(method = "lm", color = "red", se = FALSE)  # 線形回帰線を追加（任意）
#################################################################################################
head(umap_df_tcga)
# 発現量の比較
library(ggplot2)
library(reshape2)

# 悪性度関連マーカー（例: CEA, MYC, CCND1）
malignancy_markers <- c("ENSG00000170956", "ENSG00000136997", "ENSG00000110092")  # CEA, MYC, CCND1

# 発現データの取得（バージョン付きの Ensembl ID を考慮）
marker_expr_list <- lapply(malignancy_markers, function(gene) {
  gene_full <- grep(paste0("^", gene, "\\."), rownames(normalized_counts_tcga), value = TRUE)
  if (length(gene_full) > 0) {
    return(as.numeric(normalized_counts_tcga[gene_full, rownames(umap_df_tcga)]))
  } else {
    return(rep(NA, nrow(umap_df_tcga)))  # 遺伝子が見つからなかった場合
  }
})

# データフレーム化
marker_expr_df <- data.frame(
  Sample = rownames(umap_df_tcga),
  Condition = umap_df_tcga$condition,
  CDK9_ICE1_Top10 = ifelse(umap_df_tcga$CDK9_ICE1_ratio_log >= quantile(umap_df_tcga$CDK9_ICE1_ratio_log, 0.90, na.rm = TRUE), "Top10", "Other"),
  CEA = marker_expr_list[[1]],
  MYC = marker_expr_list[[2]],
  CCND1 = marker_expr_list[[3]]
)

# データをロングフォーマットに変換
marker_expr_long <- melt(marker_expr_df, id.vars = c("Sample", "Condition", "CDK9_ICE1_Top10"), 
                         variable.name = "Gene", value.name = "Expression")

# 箱ひげ図の作成
ggplot(marker_expr_long, aes(x = Gene, y = Expression, fill = CDK9_ICE1_Top10)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Expression of Malignancy Markers in CDK9/ICE1 High vs. Other",
       x = "Gene",
       y = "Expression Level",
       fill = "Group") +
  scale_fill_manual(values = c("Top10" = "red", "Other" = "gray"))
#####################################################################################################
# CEA 上位 10% の閾値を計算
threshold_cea_top10 <- quantile(umap_df_tcga$CEA_expression, 0.90, na.rm = TRUE)

# グループ分類（上位 10% を "High"、それ以外を "Other"）
umap_df_tcga$CEA_highlight <- ifelse(umap_df_tcga$CEA_expression >= threshold_cea_top10, "High", "Other")

# **赤を最前面にするために、"Other" を先に、"High" を後にプロット**
umap_df_tcga <- umap_df_tcga[order(umap_df_tcga$CEA_highlight, decreasing = TRUE), ]

umap_cea_plot <- ggplot(umap_df_tcga, aes(x = UMAP1, y = UMAP2, color = CEA_highlight)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("High" = "red", "Other" = "gray")) +
  theme_minimal() +
  labs(title = "UMAP: CEA High (Top 10%) vs. Others",
       color = "CEA Group")

# プロットの表示
print(umap_cea_plot)
ggsave("UMAP_CEA_Top10.pdf", 
       plot = umap_cea_plot,  
       width = 7, height = 5, dpi = 300, device = cairo_pdf)
#########################################################################################################]
# 癌と健常の UMAP 範囲を定義
umap_df_tcga$group <- ifelse(
  umap_df_tcga$UMAP1 >= -4 & umap_df_tcga$UMAP1 <= 3 & umap_df_tcga$UMAP2 >= -5 & umap_df_tcga$UMAP2 <= 2.5, "Tumor",
  ifelse(umap_df_tcga$UMAP1 >= 7.5 & umap_df_tcga$UMAP1 <= 10 & umap_df_tcga$UMAP2 >= 14 & umap_df_tcga$UMAP2 <= 20, "Normal", NA)
)

# NA（分類できなかったサンプル）を除外
umap_filtered <- umap_df_tcga[!is.na(umap_df_tcga$group), ]
# 癌と健常のサンプルを取得
tumor_samples <- rownames(umap_filtered[umap_filtered$group == "Tumor", ])
normal_samples <- rownames(umap_filtered[umap_filtered$group == "Normal", ])

# 正規化発現データをフィルタリング（バージョン付きの ENSEMBL ID に対応）
tumor_expr <- normalized_counts_tcga[, tumor_samples]
normal_expr <- normalized_counts_tcga[, normal_samples]

# 遺伝子ごとの平均発現値を計算
tumor_means <- rowMeans(tumor_expr, na.rm = TRUE)
normal_means <- rowMeans(normal_expr, na.rm = TRUE)

# 差分を計算
diff_expr <- data.frame(
  Gene = rownames(normalized_counts_tcga),
  Tumor_Mean = tumor_means,
  Normal_Mean = normal_means,
  Log2FC = log2((tumor_means + 1) / (normal_means + 1))  # Log2 Fold Change
)

# Fold Change の大きい順に並べる
diff_expr <- diff_expr[order(-abs(diff_expr$Log2FC)), ]
head(diff_expr, 20)  # 上位20遺伝子を表示
library(DESeq2)

# DESeq2 用データフレーム作成
coldata <- data.frame(
  row.names = c(tumor_samples, normal_samples),
  condition = factor(rep(c("Tumor", "Normal"), c(length(tumor_samples), length(normal_samples))))
)

# DESeq2 オブジェクト作成
dds <- DESeqDataSetFromMatrix(
  countData = normalized_counts_tcga[, c(tumor_samples, normal_samples)],
  colData = coldata,
  design = ~ condition
)

# 差次的発現解析の実行
dds <- DESeq(dds)

# 結果を取得
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))

# 有意な遺伝子を p 値で並べる
res_ordered <- res[order(res$pvalue), ]
head(res_ordered, 20)  # 上位20遺伝子を表示
library(ggplot2)

# 結果をデータフレームに変換
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)
res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 3, "Significant", "Not Significant")

# Volcano Plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "Log2 Fold Change",
       y = "-log10(p-value)")
getwd()
# CSV ファイルとして保存
write.csv(diff_expr, file = "diff_expr_top20.csv", row.names = FALSE)

# res_ordered も保存
write.csv(res_ordered, file = "res_ordered_top20.csv", row.names = FALSE)
library(ggplot2)
library(ggrepel)  # ラベル表示用

# ICE1 と CDK9 の遺伝子 ID
ice1_gene <- "ENSG00000164151"
cdk9_gene <- "ENSG00000136807"

# ICE1/CDK9 の発現データを抽出
highlight_genes <- res_df[res_df$Gene %in% c(ice1_gene, cdk9_gene), ]

res_df$Gene <- gsub("\\.[0-9]+$", "", res_df$Gene)
# NA を含む行を削除
res_df_clean <- na.omit(res_df)
# NA を含む行をすべて削除
res_df_clean <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$pvalue) & !is.na(res_df$padj), ]

# Volcano Plot
volcano_plot <- ggplot(res_df_clean, aes(x = log2FoldChange, y = -log10(pvalue), color = Significant)) +
  geom_point(alpha = 0.6) +  # すべての点を描画
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Gene Expression",
       x = "Log2 Fold Change",
       y = "-log10(p-value)") +
  
  # ICE1 と CDK9 を青で最前面に描画
  geom_point(data = highlight_genes, aes(x = log2FoldChange, y = -log10(pvalue)), 
             color = "blue", size = 6, shape = 17) +  # ICE1/CDK9 を大きく
  geom_text_repel(data = highlight_genes, aes(x = log2FoldChange, y = -log10(pvalue), label = Gene), 
                  size = 5, color = "blue", fontface = "bold")  # ラベルを表示

# プロットを表示
print(volcano_plot)
print(ice1_gene %in% res_df$Gene)  # TRUE なら OK
print(cdk9_gene %in% res_df$Gene)  # TRUE なら OK
print(sum(is.na(res_df$log2FoldChange)))  # log2FoldChange の NA の数
print(sum(is.na(res_df$pvalue)))         # pvalue の NA の数
res_df_clean <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$pvalue), ]
res_df_clean <- res_df_clean[!is.na(res_df_clean$padj), ]
print(highlight_genes)
##############################################################################################################
# CDK9/ICE1 の比率データを取得
cdk9_ice1_ratio <- umap_filtered$CDK9_ICE1_ratio_log

# 上位10%のカットオフ値を計算
cutoff_value_1 <- quantile(cdk9_ice1_ratio, 0.9)

# カットオフ値を表示
print(cutoff_value_1)
# 癌 (Tumor) と正常 (Normal) のラベル
true_labels <- umap_filtered$group  # 実際の分類 (Tumor or Normal)

# CDK9/ICE1 の比率に基づく予測ラベル
predicted_labels_1 <- ifelse(umap_filtered$CDK9_ICE1_ratio_log >= cutoff_value_1, "Tumor", "Normal")

# 混同行列を作成
conf_matrix_1 <- table(Actual = true_labels, Predicted = predicted_labels_1)

# 混同行列を表示
print(conf_matrix_1)

library(ggplot2)
library(gridExtra)

# 混同行列の各要素
TP_1 <- conf_matrix_1["Tumor", "Tumor"]   # 真陽性
FN_1 <- conf_matrix_1["Tumor", "Normal"]  # 偽陰性
FP_1 <- conf_matrix_1["Normal", "Tumor"]  # 偽陽性
TN_1 <- conf_matrix_1["Normal", "Normal"] # 真陰性

# 指標の計算
sensitivity_1 <- TP_1 / (TP_1 + FN_1)  # 感度
specificity_1 <- TN_1 / (TN_1 + FP_1)  # 特異度
FPR <- FP / (FP + TN)          # 偽陽性率
FNR <- FN / (TP + FN)          # 偽陰性率

# ヒートマップ用データフレーム
conf_matrix_df_1 <- as.data.frame(as.table(conf_matrix_1))

# ✅ **正方形のヒートマップ**
heatmap_plot <- ggplot(conf_matrix_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 4, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  coord_fixed(ratio = 1) +  # ✅ **正方形にする**
  labs(title = "CDK9/ICE1 matrix")

# ✅ **4つの指標を含む棒グラフ**
metrics_text <- data.frame(
  Metric = c("Sensitivity", "Specificity"),
  Value = c(sensitivity_1, specificity_1)
)

metrics_plot <- ggplot(metrics_text, aes(x = Metric, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Sensitivity" = "blue", 
                               "Specificity" = "green")) +
  theme_minimal() +
  labs(title = "CDK9/ICE1", y = "Value") +
  geom_text(aes(label = round(Value, 3)), vjust = -0.5, size = 3, fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0,1)
# **X軸のラベルを斜めにする**

# ✅ **2つの図を統合**
grid.arrange(heatmap_plot, metrics_plot, ncol = 2)

# ✅ CEA のカットオフ値を計算
cutoff_cea <- quantile(umap_filtered$CEA_expression, 0.9)

# ✅ 実際のラベル (True labels)
true_labels <- umap_filtered$group  # Tumor or Normal

# ✅ CEA の発現量に基づく予測ラベル
predicted_labels <- ifelse(umap_filtered$CEA_expression >= cutoff_cea, "Tumor", "Normal")

# ✅ 混同行列を作成
conf_matrix_cea <- table(Actual = true_labels, Predicted = predicted_labels)

# ✅ 混同行列の各要素を取得
TP_1 <- conf_matrix_cea["Tumor", "Tumor"]   # 真陽性
FN_1 <- conf_matrix_cea["Tumor", "Normal"]  # 偽陰性
FP_1 <- conf_matrix_cea["Normal", "Tumor"]  # 偽陽性
TN_1 <- conf_matrix_cea["Normal", "Normal"] # 真陰性

# ✅ 感度・特異度・偽陽性率・偽陰性率を計算
sensitivity_1 <- TP / (TP + FN)  # 感度 (Sensitivity)
specificity_1 <- TN / (TN + FP)  # 特異度 (Specificity)
FPR_1 <- FP / (FP + TN)          # 偽陽性率 (False Positive Rate)
FNR_1 <- FN / (TP + FN)          # 偽陰性率 (False Negative Rate)

# ✅ ヒートマップ用データフレーム
conf_matrix_df <- as.data.frame(as.table(conf_matrix_cea))

# ✅ **正方形ヒートマップ**
heatmap_plot <- ggplot(conf_matrix_df, aes(x = Predicted, y = Actual, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 4, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  coord_fixed(ratio = 1) +  # ✅ **正方形にする**
  labs(title = "CEA")

# ✅ **4つの指標を含む棒グラフ**
metrics_text <- data.frame(
  Metric = c("Sensitivity", "Specificity"),
  Value = c(sensitivity_1, specificity_1)
)

metrics_plot <- ggplot(metrics_text, aes(x = Metric, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Sensitivity" = "blue", 
                               "Specificity" = "green" 
                               )) +
  theme_minimal() +
  labs(title = "CEA", y = "Value") +
  geom_text(aes(label = round(Value, 3)), vjust = -0.5, size = 3, fontface = "bold") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0,1)  # **X軸のラベルを斜めにする**

# ✅ **2つの図を統合**
final_plot <- grid.arrange(heatmap_plot, metrics_plot, ncol = 2)

# ✅ **PDFで保存 (7:5 の比率)**
ggsave("CEA_Confusion_Matrix_and_Metrics.pdf", plot = final_plot, width = 7, height = 5)

install.packages("pROC")
library(pROC)

# ✅ 実際のラベル (Ground Truth)
true_labels <- ifelse(umap_filtered$group == "Tumor", 1, 0)  # 1: Tumor, 0: Normal

# ✅ 予測スコア (CDK9/ICE1 の比)
predicted_scores <- umap_filtered$CDK9_ICE1_ratio

# ✅ ROC 曲線を計算
roc_curve <- roc(true_labels, predicted_scores)

# ✅ AUC を取得
auc_value <- auc(roc_curve)

# ✅ ROC 曲線を ggplot2 で描画
roc_df <- data.frame(TPR = rev(roc_curve$sensitivities), 
                     FPR = rev(1 - roc_curve$specificities))

roc_plot <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "blue", size = 1) +
  geom_abline(linetype = "dashed", color = "gray") +
  theme_minimal() +
  labs(title = "ROC Curve for CDK9/ICE1",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  annotate("text", x = 0.7, y = 0.2, label = paste("AUC =", round(auc_value, 3)), 
           color = "red", size = 6, fontface = "bold")

# ✅ PDFで保存 (7:5 の比率)
ggsave("CDK9_ICE1_ROC_Curve.pdf", plot = roc_plot, width = 7, height = 5)
# ✅ ROC曲線の計算
roc_curve <- roc(umap_filtered$group, umap_filtered$CDK9_ICE1_ratio, levels = c("Normal", "Tumor"))

# ✅ Youden’s Index を使って最適なしきい値を計算
roc_data <- data.frame(
  Sensitivity = roc_curve$sensitivities,
  Specificity = roc_curve$specificities,
  Threshold = roc_curve$thresholds
)
roc_data$Youden_Index <- roc_data$Sensitivity + roc_data$Specificity - 1

# ✅ 最大の Youden’s Index を持つしきい値を取得
optimal_threshold <- roc_data$Threshold[which.max(roc_data$Youden_Index)]
print(optimal_threshold)
##################################################################################################################
#CDK9/ICE1のｶｯﾄｵﾌをした時の特徴の違い
# ✅ CDK9/ICE1 の比のカットオフ (上位10%の閾値を求める)
cutoff_value_1 <- quantile(umap_filtered$CDK9_ICE1_ratio, 0.9)

# ✅ カットオフを基準に "High" と "Low" に分類
umap_filtered$CDK9_ICE1_group <- ifelse(umap_filtered$CDK9_ICE1_ratio >= cutoff_value_1, "High", "Low")

# ✅ 高い群・低い群のサンプルを取得
high_group_samples <- rownames(umap_filtered[umap_filtered$CDK9_ICE1_group == "High", ])
low_group_samples <- rownames(umap_filtered[umap_filtered$CDK9_ICE1_group == "Low", ])

# ✅ グループごとのサンプル数を確認
table(umap_filtered$CDK9_ICE1_group)
library(DESeq2)

# ✅ DESeq2 のデータセット作成
coldata <- data.frame(
  row.names = c(high_group_samples, low_group_samples),
  group = factor(rep(c("High", "Low"), c(length(high_group_samples), length(low_group_samples))))
)
# 数値データを整数に丸める
normalized_counts_tcga_int <- round(normalized_counts_tcga)

# 再度DESeq2用データセットを作成
dds <- DESeqDataSetFromMatrix(
  countData = normalized_counts_tcga_int[, c(high_group_samples, low_group_samples)],
  colData = coldata,
  design = ~ group
)

# ✅ 差次的発現解析の実行
dds <- DESeq(dds)

# ✅ 結果を取得
res <- results(dds, contrast = c("group", "High", "Low"))

# ✅ 有意な遺伝子を p 値で並べる
res_ordered <- res[order(res$pvalue), ]
head(res_ordered, 20)  # 上位20遺伝子を表示

# ✅ 結果を CSV に保存
write.csv(as.data.frame(res_ordered), "CDK9_ICE1_High_vs_Low_DEGs.csv")
library(ggrepel)

# ✅ ICE1, CDK9 の遺伝子 ID
ice1_gene <- "ENSG00000164151"
cdk9_gene <- "ENSG00000136807"

# ✅ データフレームに変換
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)
res_df$Significant <- ifelse(res_df$padj < 0.0005 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")

# ✅ ICE1, CDK9 を青色で強調表示
res_df$Highlight <- ifelse(res_df$Gene %in% c(ice1_gene, cdk9_gene), "Highlight", res_df$Significant)

# ✅ Volcano Plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Highlight)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray", "Highlight" = "blue")) +
  geom_text_repel(data = res_df[res_df$Gene %in% c(ice1_gene, cdk9_gene), ],
                  aes(label = ifelse(Gene == ice1_gene, "ICE1", "CDK9")),
                  size = 5, color = "blue", fontface = "bold") +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Gene Expression (CDK9/ICE1 High vs Low)",
       x = "Log2 Fold Change",
       y = "-log10(p-value)")

# ✅ PDF に保存 (7:5 の比率)
ggsave("CDK9_ICE1_VolcanoPlot_with_Labels.pdf", plot = volcano_plot, width = 7, height = 5)

# ✅ 有意に上昇した遺伝子を抽出
upregulated_genes <- res_df$Gene[res_df$log2FoldChange > 1 & res_df$padj < 0.005]

# ✅ 結果を CSV に保存
write.csv(upregulated_genes, "Upregulated_Genes_CDK9_ICE1_High.csv", row.names = FALSE)

# ✅ 上昇遺伝子の数を確認
length(upregulated_genes)
library(clusterProfiler)
library(org.Hs.eg.db)  # ヒトの遺伝子アノテーション

# ✅ ENSG ID を ENTREZ ID に変換
upregulated_entrez <- bitr(upregulated_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# ✅ GO Enrichment 解析 (Biological Process)
go_enrich <- enrichGO(
  gene = upregulated_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)



# ✅ 結果を CSV に保存
write.csv(as.data.frame(go_enrich), "GO_Enrichment_CDK9_ICE1_High.csv", row.names = FALSE)

# ✅ バブルプロットを作成
go_plot <- dotplot(go_enrich, showCategory = 10, title = "GO Enrichment (CDK9/ICE1 High Group)")

# ✅ PDF に保存 (7:5 の比率)
ggsave("GO_Enrichment_BubblePlot_CDK9_ICE1_High.pdf", plot = go_plot, width = 7, height = 5)

library(clusterProfiler)
library(org.Hs.eg.db)  # ヒトの遺伝子アノテーション

# ✅ 減少している遺伝子（Log2FC < -1 & padj < 0.005）を抽出
downregulated_genes_cdk9 <- res_df$Gene[res_df$log2FoldChange < -1 & res_df$padj < 0.005]

# ✅ 結果を CSV に保存
write.csv(downregulated_genes_cdk9, "Downregulated_Genes_CDK9_ICE1_Low.csv", row.names = FALSE)

# ✅ 減少遺伝子の数を確認
length(downregulated_genes_cdk9)

# ✅ ENSG ID を ENTREZ ID に変換
downregulated_entrez_cdk9 <- bitr(downregulated_genes_cdk9, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# ✅ GO Enrichment 解析 (Biological Process)
go_enrich_down_cdk9 <- enrichGO(
  gene = downregulated_entrez_cdk9$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# ✅ 結果を CSV に保存
write.csv(as.data.frame(go_enrich_down_cdk9), "GO_Enrichment_CDK9_ICE1_Low.csv", row.names = FALSE)

# ✅ バブルプロットを作成
go_plot_down_cdk9 <- dotplot(go_enrich_down_cdk9, showCategory = 10, title = "GO Enrichment (CDK9/ICE1 Low Group)")

# ✅ PDF に保存 (7:5 の比率)
ggsave("GO_Enrichment_BubblePlot_CDK9_ICE1_Low.pdf", plot = go_plot_down_cdk9, width = 7, height = 5)


# ✅ DESeq2 の結果から上位5遺伝子を選択
top12_genes <- rownames(res_ordered)[1:12]

# ✅ 上位512伝子の発現をプロット
gene_expr_long <- reshape2::melt(normalized_counts_tcga[top12_genes, c(high_group_samples, low_group_samples)])
colnames(gene_expr_long) <- c("Gene", "Sample", "Expression")
gene_expr_long$Group <- umap_filtered$CDK9_ICE1_group[match(gene_expr_long$Sample, rownames(umap_filtered))]

pICE <- ggplot(gene_expr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin() +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_minimal() +
  labs(title = "Top 12 Differentially Expressed Genes by CDK9/ICE1 Group", x = "Group", y = "Expression")
ggsave("CDK9_ICE1_ViolinPlot_with_Labels.pdf", plot = pICE, width = 7, height = 5)
write.csv(as.data.frame(top12_genes), "CDK9_ICE1_top12genes.csv", row.names = FALSE)
###########################################################################################################
# ✅ CEA のカットオフ (上位10%の閾値)
cutoff_value_cea <- quantile(umap_filtered$CEA_expression, 0.9)

# ✅ CEA の発現量を基準に「High」と「Low」に分類
umap_filtered$CEA_group <- ifelse(umap_filtered$CEA_expression >= cutoff_value_cea, "High", "Low")

# ✅ 高群・低群のサンプルを取得
high_group_cea_samples <- rownames(umap_filtered[umap_filtered$CEA_group == "High", ])
low_group_cea_samples <- rownames(umap_filtered[umap_filtered$CEA_group == "Low", ])

# ✅ グループごとのサンプル数を確認
table(umap_filtered$CEA_group)

library(DESeq2)

# ✅ CEA High/Low 用の colData 作成
coldata_cea <- data.frame(
  row.names = c(high_group_cea_samples, low_group_cea_samples),
  group = factor(rep(c("High", "Low"), c(length(high_group_cea_samples), length(low_group_cea_samples))))
)

# 数値データを整数に丸める
normalized_counts_tcga_int <- round(normalized_counts_tcga)

# DESeq2データセット作成
dds_cea <- DESeqDataSetFromMatrix(
  countData = normalized_counts_tcga_int[, c(high_group_cea_samples, low_group_cea_samples)],
  colData = coldata_cea,
  design = ~ group
)

# 差次的発現解析の実行
dds_cea <- DESeq(dds_cea)

# 結果を取得
res_cea <- results(dds_cea, contrast = c("group", "High", "Low"))

# 有意な遺伝子を p 値で並べる
res_ordered_cea <- res_cea[order(res_cea$pvalue), ]
head(res_ordered_cea, 20)  # 上位20遺伝子を表示

# 結果を CSV に保存
write.csv(as.data.frame(res_ordered_cea), "CEA_High_vs_Low_DEGs.csv")

library(ggplot2)
library(ggrepel)

# CEAの遺伝子 ID
cea_gene <- "ENSG00000170956"

# ✅ Volcano Plot データフレームを作成
res_df_cea <- as.data.frame(res_cea)
res_df_cea$Gene <- rownames(res_df_cea)
res_df_cea$Significant <- ifelse(res_df_cea$padj < 0.0005 & abs(res_df_cea$log2FoldChange) > 1, "Significant", "Not Significant")

# ✅ CEA を青色で強調表示
res_df_cea$Highlight <- ifelse(res_df_cea$Gene == cea_gene, "Highlight", res_df_cea$Significant)

# ✅ Volcano Plot
volcano_plot_cea <- ggplot(res_df_cea, aes(x = log2FoldChange, y = -log10(pvalue), color = Highlight)) +
  geom_point(alpha = 1) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray", "Highlight" = "blue")) +
  geom_text_repel(data = res_df_cea[res_df_cea$Gene == cea_gene, ],
                  aes(label = "CEA"),
                  size = 5, color = "blue", fontface = "bold") +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Gene Expression (CEA High vs Low)",
       x = "Log2 Fold Change",
       y = "-log10(p-value)")

# ✅ PDF に保存 (7:5 の比率)
ggsave("CEA_VolcanoPlot_with_Labels.pdf", plot = volcano_plot_cea, width = 7, height = 5)

library(clusterProfiler)
library(org.Hs.eg.db)

# 有意に上昇した遺伝子を抽出
upregulated_genes_cea <- res_df_cea$Gene[res_df_cea$log2FoldChange > 1 & res_df_cea$padj < 0.005]

# ENSG ID を ENTREZ ID に変換
upregulated_entrez_cea <- bitr(upregulated_genes_cea, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO Enrichment 解析 (Biological Process)
go_enrich_cea <- enrichGO(
  gene = upregulated_entrez_cea$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# 結果を CSV に保存
write.csv(as.data.frame(go_enrich_cea), "GO_Enrichment_CEA_High.csv", row.names = FALSE)

# GO Enrichment のバブルプロット
go_plot_cea <- dotplot(go_enrich_cea, showCategory = 10, title = "GO Enrichment (CEA High Group)")

# PDF に保存 (7:5 の比率)
ggsave("GO_Enrichment_BubblePlot_CEA_High.pdf", plot = go_plot_cea, width = 7, height = 5)

library(clusterProfiler)
library(org.Hs.eg.db)  # ヒトの遺伝子アノテーション

# ✅ 減少している遺伝子（Log2FC < -1 & padj < 0.005）を抽出
downregulated_genes <- res_df_cea$Gene[res_df_cea$log2FoldChange < -1 & res_df_cea$padj < 0.005]

# ✅ 結果を CSV に保存
write.csv(downregulated_genes, "Downregulated_Genes_CEA_Low.csv", row.names = FALSE)

# ✅ 減少遺伝子の数を確認
length(downregulated_genes)

# ✅ ENSG ID を ENTREZ ID に変換
downregulated_entrez <- bitr(downregulated_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# ✅ GO Enrichment 解析 (Biological Process)
go_enrich_down <- enrichGO(
  gene = downregulated_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# ✅ 結果を CSV に保存
write.csv(as.data.frame(go_enrich_down), "GO_Enrichment_CEA_Low.csv", row.names = FALSE)

# ✅ バブルプロットを作成
go_plot_down <- dotplot(go_enrich_down, showCategory = 10, title = "GO Enrichment (CEA Low Group)")

# ✅ PDF に保存 (7:5 の比率)
ggsave("GO_Enrichment_BubblePlot_CEA_Low.pdf", plot = go_plot_down, width = 7, height = 5)




##############################################################################################################
library(reshape2)

# 上位12遺伝子を抽出
top12_genes_cea <- rownames(res_ordered_cea)[1:12]

# データをロングフォーマットに変換
gene_expr_long_cea <- melt(normalized_counts_tcga[top12_genes_cea, c(high_group_cea_samples, low_group_cea_samples)])
colnames(gene_expr_long_cea) <- c("Gene", "Sample", "Expression")
gene_expr_long_cea$Group <- umap_filtered$CEA_group[match(gene_expr_long_cea$Sample, rownames(umap_filtered))]

# バイオリンプロット
pCEA <- ggplot(gene_expr_long_cea, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin() +
  facet_wrap(~ Gene, scales = "free_y") +
  theme_minimal() +
  labs(title = "Top 12 Differentially Expressed Genes by CEA Group", x = "Group", y = "Expression")

# PDF に保存 (7:5 の比率)
ggsave("CEA_ViolinPlot_with_Labels.pdf", plot = pCEA, width = 7, height = 5)

# 上位12遺伝子リストを CSV に保存
write.csv(as.data.frame(top12_genes_cea), "CEA_top12genes.csv", row.names = FALSE)
