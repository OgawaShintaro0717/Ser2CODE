"C:/Documents/20241117RNAseqIntergration/20221021_CoilinChAPpeakCall_Annotated_1.csv"
data <- read.csv("C:/Documents/20241117RNAseqIntergration/expented_count.csv")
head(data)
library(biomaRt)
# BiomaRtに接続
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# ENSG IDに基づいて遺伝子情報を取得
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = data$gene_id,
  mart = mart
)

# 元のデータフレームに情報を結合
data_with_gene_info <- merge(
  data, gene_info,
  by.x = "gene_id", by.y = "ensembl_gene_id",
  all.x = TRUE
)

# 確認
head(data_with_gene_info)
#####################################################################################
CB <- read.csv("C:/Documents/20241117RNAseqIntergration/20221021_CoilinChAPpeakCall_Annotated_1.csv")
head(CB)
# 必要な列を抽出
data_1_ENSG <- unique(data_1$feature)

# 新しい列を追加
data$CB <- data$gene_id %in% data_1_ENSG

# 結果を確認
head(data)
subset(data, gene_id == "ENSG00000206652")
#####################################################################################
#DFの作り直し
# 1, 16, 17, 18, 2~8列を抽出してDF_7SKを作成
DF_7SK <- data_with_gene_info[, c(1, 16, 17, 18, 2:9)]

# 1, 16, 17, 18, 9~16列を抽出してDF_AFF4を作成
DF_AFF4 <- data_with_gene_info[, c(1, 16, 17, 18, 10:16)]

# 作成したデータフレームの確認
head(DF_7SK)
head(DF_AFF4)
write.csv(DF_7SK, file = "C:/Documents/20241117RNAseqIntergration/DF_7SK.csv", row.names = FALSE)
write.csv(DF_AFF4, file = "C:/Documents/20241117RNAseqIntergration/DF_AFF4.csv", row.names = FALSE)

DF_7SK <- read.csv("C:/Documents/20241117RNAseqIntergration/DF_7SK.csv")
DF_AFF4 <- read.csv("C:/Documents/20241117RNAseqIntergration/DF_AFF4.csv")
head(DF_7SK)
head(DF_AFF4)

###################################################################################
DF_AFF4_p05 <- subset(DF_AFF4, p_value < 0.05 & is.finite(log2FC))
head(DF_AFF4_p05)
# log2FCの絶対値を計算して新しい列を作成
DF_AFF4_p05$abs_log2FC <- abs(DF_AFF4_p05$log2FC)

# log2FCの絶対値が大きい順にソート
DF_AFF4_p05_sorted <- DF_AFF4_p05[order(-DF_AFF4_p05$abs_log2FC), ]
# 上位500個の遺伝子を取得
DF_AFF4_top500 <- head(DF_AFF4_p05_sorted, 1000)
# 結果を表示
head(DF_AFF4_top500)

# gene_biotypeの頻度を計算
biotype_counts <- table(DF_AFF4_top500$gene_biotype)

# 全体に占める割合を計算
biotype_percentages <- prop.table(biotype_counts) * 100

# biotype_countsをデータフレームに変換
biotype_df <- data.frame(
  Category = names(biotype_counts),
  Count = as.numeric(biotype_counts)
)

threshold <- 2
# 2%以下の項目をまとめるための分類
biotype_df$Category <- ifelse(
  biotype_percentages > threshold, 
  biotype_df$Category, 
  "Others"
)

# "Others"を含むカテゴリごとに集計
biotype_adjusted_counts <- aggregate(
  Count ~ Category, 
  data = biotype_df, 
  FUN = sum
)

# 全体に占める割合を再計算
biotype_adjusted_percentages <- prop.table(biotype_adjusted_counts$Count) * 100
# 色を青系に設定
blue_shades <- colorRampPalette(c("lightblue", "blue"))(length(biotype_adjusted_percentages))

# 円グラフを描画
pie(biotype_adjusted_percentages,
    labels = paste0(biotype_adjusted_counts$Category, " (", round(biotype_adjusted_percentages, 1), "%)"),
    main = "Gene Biotype Distribution",
    col = blue_shades,
    cex = 0.8) # ラベルの文字サイズを調整
###################################################################################
DF_7SK_p05 <- subset(DF_7SK, p_value < 0.05 & is.finite(log2FC))
head(DF_7SK_p05)
# log2FCの絶対値を計算して新しい列を作成
DF_7SK_p05$abs_log2FC <- abs(DF_7SK_p05$log2FC)

# log2FCの絶対値が大きい順にソート
DF_7SK_p05_sorted <- DF_7SK_p05[order(-DF_7SK_p05$abs_log2FC), ]
# 上位500個の遺伝子を取得
DF_7SK_top500 <- head(DF_7SK_p05_sorted, 1000)
# 結果を表示
head(DF_7SK_top500)

# gene_biotypeの頻度を計算
biotype_counts <- table(DF_7SK_top500$gene_biotype)

# 全体に占める割合を計算
biotype_percentages <- prop.table(biotype_counts) * 100

# biotype_countsをデータフレームに変換
biotype_df <- data.frame(
  Category = names(biotype_counts),
  Count = as.numeric(biotype_counts)
)

threshold <- 2
# 2%以下の項目をまとめるための分類
biotype_df$Category <- ifelse(
  biotype_percentages > threshold, 
  biotype_df$Category, 
  "Others"
)

# "Others"を含むカテゴリごとに集計
biotype_adjusted_counts <- aggregate(
  Count ~ Category, 
  data = biotype_df, 
  FUN = sum
)

# 全体に占める割合を再計算
biotype_adjusted_percentages <- prop.table(biotype_adjusted_counts$Count) * 100
# 色を青系に設定
blue_shades <- colorRampPalette(c("lightblue", "blue"))(length(biotype_adjusted_percentages))

# 円グラフを描画
pie(biotype_adjusted_percentages,
    labels = paste0(biotype_adjusted_counts$Category, " (", round(biotype_adjusted_percentages, 1), "%)"),
    main = "Gene Biotype Distribution",
    col = blue_shades,
    cex = 0.8) # ラベルの文字サイズを調整

###################################################################################
DF_7SK <- read.csv("C:/Documents/20241117RNAseqIntergration/DF_7SK.csv")
head(DF_7SK)

# 1. siCTとsi7SKの平均を計算
DF_7SK$siCT_mean <- rowMeans(DF_7SK[, 5:8])  # sample1〜sample4
DF_7SK$si7SK_mean <- rowMeans(DF_7SK[, 9:12])  # sample5〜sample8
DF_7SK[, 5:12] <- lapply(DF_7SK[, 5:12], as.numeric)
# 2. Fold Change (FC) を計算
DF_7SK$log2FC <- log2(DF_7SK$si7SK_mean / DF_7SK$siCT_mean)

# 3. p値の計算 (t検定)
p_values <- apply(DF_7SK[, 5:12], 1, function(x) {
  tryCatch({
    t.test(x[1:4], x[5:8])$p.value
  }, error = function(e) {
    return(NA)  # エラーが発生した場合はNAを返す
  })
})
DF_7SK$p_value <- p_values

# 4. -log10(p-value)を計算して、Volcano plotに使用する
DF_7SK$log_p_value <- -log10(DF_7SK$p_value)

# NAやInfを含む行を削除
DF_cleaned <- DF_7SK[!is.na(DF_7SK$log2FC) & !is.infinite(DF_7SK$log2FC) &
                       !is.na(DF_7SK$p_value), ]
head(DF_cleaned)
# 5. Volcano plotの作成
library(ggplot2)

# CBがTRUEの点を赤、FALSEの点を灰色で表示
volcano <- ggplot(DF_cleaned, aes(x=log2FC, y=-log10(p_value))) +
  # 灰色の点を先に描画
  geom_point(data = subset(DF_cleaned, CB == FALSE), aes(color = CB), size = 2) +
  # 赤い点を後に描画
  geom_point(data = subset(DF_cleaned, CB == TRUE), aes(color = CB), size = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +  # 色の指定
  theme_minimal() +
  labs(x="log2 Fold Change", y="-log10(p-value)") +
  xlim(-5, 5) +   
  ylim(0, 7) +   
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
ggsave("C:/Documents/20241117RNAseqIntergration/si7SK_bigwigVolcano.pdf", plot = volcano, width = 10, height = 8)
# CBがTRUEの点の数を数える
sum(DF_cleaned$CB == TRUE)
#################################################################################################################
# 必要なライブラリをロード
library(ggplot2)

# p-valueが0.05未満であるか、log2FCの絶対値が上位1000位以内のものをフィルタリング
DF_cleaned$log_p_value <- -log10(DF_cleaned$p_value)
DF_cleaned$abs_log2FC <- abs(DF_cleaned$log2FC)

# log2FCの絶対値で上位1000位を取得 (全行数が1000未満の場合はすべて含まれます)
top_1000_abs_log2FC <- sort(DF_cleaned$abs_log2FC, decreasing = TRUE)[min(1000, nrow(DF_cleaned))]

# 条件を満たすかどうかを判定
DF_cleaned$highlight <- ifelse(DF_cleaned$p_value < 0.05 & DF_cleaned$abs_log2FC >= top_1000_abs_log2FC, "Significant", "Not Significant")

# 遺伝子名がRN7SKのデータポイントを特定してハイライト
DF_cleaned$color <- ifelse(DF_cleaned$external_gene_name == "RN7SK", "RN7SK", DF_cleaned$highlight)

# Volcano Plotを描画
volcano_plot <- ggplot(DF_cleaned, aes(x = log2FC, y = log_p_value, color = color)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("Not Significant" = "black", "Significant" = "red", "RN7SK" = "blue")) +
  labs(
    title = "si7SK RNAseq p_value < 0.05 & top_1000_abs_log2FC genes",
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  ) +
  xlim(-6, 6) +   
  ylim(0, 7) +   
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )
# PDFに保存
output_path <- "C:/Documents/20241117RNAseqIntergration/volcano_plot_RN7SK_KD_p_value005&top_1000_abs_log2FCgenes.pdf"
pdf(output_path, width = 8, height = 6) # サイズは適宜調整
print(volcano_plot)
dev.off()

# 保存完了メッセージ
cat("PDFファイルが以下に保存されました:", output_path, "\n")
#################################################################################################################################

# RN7SKのデータポイントを抽出
RN7SK_row <- DF_cleaned[DF_cleaned$external_gene_name == "RN7SK", ]
other_rows <- DF_cleaned[DF_cleaned$external_gene_name != "RN7SK", ]

# Volcano Plotを描画
volcano_plot_7SK <- ggplot() +
  # 他の点を描画 (最初に描画する)
  geom_point(data = other_rows, aes(x = log2FC, y = log_p_value, color = color), alpha = 1) +
  # RN7SKの点を描画 (最前面に描画)
  geom_point(data = RN7SK_row, aes(x = log2FC, y = log_p_value), color = "blue", size = 3, alpha = 1) +
  # RN7SKにラベルを追加 (左上に配置)
  geom_text(data = RN7SK_row, aes(x = log2FC, y = log_p_value, label = external_gene_name), 
            hjust = 1.2, vjust = -0.2, size = 5, color = "blue") +
  scale_color_manual(values = c("Not Significant" = "black", "Significant" = "red")) +
  labs(
    title = "7SK RNAseq p_value < 0.05 & top_1000_abs_log2FC genes",
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  ) +
  xlim(-6, 6) +   
  ylim(0, 7) +   
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

# PDFに保存
output_path <- "C:/Documents/20241117RNAseqIntergration/volcano_plot_7SK_KD_p_value005&top_1000_abs_log2FCgenes_with_label.pdf"
pdf(output_path, width = 8, height = 6) # サイズは適宜調整
print(volcano_plot_7SK)
dev.off()

# 保存完了メッセージ
cat("PDFファイルが以下に保存されました:", output_path, "\n")

#################################################################################################################################
# 遺伝子名がRN7SK、PANK4、AKAP12のデータポイントを特定してハイライト
DF_cleaned$color <- ifelse(
  DF_cleaned$external_gene_name == "RN7SK", "RN7SK",
  ifelse(DF_cleaned$external_gene_name == "PANK4", "PANK4",
         ifelse(DF_cleaned$external_gene_name == "AKAP12", "AKAP12",
                DF_cleaned$highlight))
)

# 点のサイズ列をカテゴリ化してPANK4とAKAP12を目立たせる
DF_cleaned$size <- factor(ifelse(
  DF_cleaned$external_gene_name %in% c("PANK4", "AKAP12"), "Large", "Small"
))

# プロットを作成
volcano_plot <- ggplot(DF_cleaned, aes(x = log2FC, y = log_p_value, color = color, size = size)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(
    "Not Significant" = "black", 
    "Significant" = "red", 
    "RN7SK" = "blue", 
    "PANK4" = "green", 
    "AKAP12" = "orange"
  )) +
  scale_size_manual(values = c("Small" = 1, "Large" = 3)) + # サイズを指定
  labs(
    title = "si7SK RNAseq p_value < 0.05 & top_1000_abs_log2FC genes",
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  ) +
  xlim(-5, 5) +   
  ylim(0, 7) +   
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

# PDFに保存
output_path <- "C:/Documents/20241117RNAseqIntergration/volcano_plot.pdf"
pdf(output_path, width = 8, height = 6) # サイズは適宜調整
print(volcano_plot)
dev.off()

# 保存完了メッセージ
cat("PDFファイルが以下に保存されました:", output_path, "\n")
##########################################################################################################
# CBがTRUEのデータのみを抽出してDF_7SK_CBに保存
DF_7SK_CB <- DF_cleaned[DF_cleaned$CB == TRUE, ]
head(DF_7SK_CB)
# log2FCとp_valueの計算
DF_7SK_CB$log2FC <- log2(DF_7SK_CB$si7SK_mean / DF_7SK_CB$siCT_mean)

# p値の計算 (t検定)
p_values <- apply(DF_7SK_CB[, 5:12], 1, function(x) {
  tryCatch({
    t.test(x[1:4], x[5:8])$p.value
  }, error = function(e) {
    return(NA)  # エラーが発生した場合はNAを返す
  })
})
DF_7SK_CB$p_value <- p_values

# -log10(p-value)を計算して、Volcano plotに使用する
DF_7SK_CB$log_p_value <- -log10(DF_7SK_CB$p_value)

# NAやInfを含む行を削除
DF_cleaned_CB <- DF_7SK_CB[!is.na(DF_7SK_CB$log2FC) & !is.infinite(DF_7SK_CB$log2FC) &
                             !is.na(DF_7SK_CB$p_value), ]

# 5. Volcano plotの作成
library(ggplot2)

# CBがTRUEの点を赤、FALSEの点を灰色で表示
volcano <- ggplot(DF_cleaned_CB, aes(x=log2FC, y=log_p_value)) +
  geom_point(aes(color=CB), size=2) +
  scale_color_manual(values = c("TRUE" = "red")) +  # CBがTRUEのものだけ赤
  theme_minimal() +
  labs(x="log2 Fold Change", y="-log10(p-value)") +
  xlim(-5, 5) +   
  ylim(0, 7) +   
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))

# Volcano plotを表示
print(volcano)

# 画像として保存
ggsave("C:/Documents/20241117RNAseqIntergration/si7SK_CB_Volcano.pdf", plot = volcano, width = 10, height = 8)



###################################################################################
DF_AFF4 <- read.csv("C:/Documents/20241117RNAseqIntergration/DF_AFF4.csv")
head(DF_AFF4)
# 1. ICE1とAFF4の平均を計算
DF_AFF4$ICE1_mean <- rowMeans(DF_AFF4[, 5:7])  # sample9〜sample11
DF_AFF4$AFF4_mean <- rowMeans(DF_AFF4[, 8:10])  # sample12〜sample13
DF_AFF4[, 5:10] <- lapply(DF_AFF4[, 5:10], as.numeric)

# 2. Fold Change (FC) を計算
DF_AFF4$log2FC <- log2(DF_AFF4$AFF4_mean / DF_AFF4$ICE1_mean)

# 3. p値の計算 (t検定)
p_values <- apply(DF_AFF4[, 5:10], 1, function(x) {
  tryCatch({
    t.test(x[1:3], x[4:6])$p.value  # ICE1とAFF4でt検定
  }, error = function(e) {
    return(NA)  # エラーが発生した場合はNAを返す
  })
})
DF_AFF4$p_value <- p_values

# 4. -log10(p-value)を計算して、Volcano plotに使用する
DF_AFF4$log_p_value <- -log10(DF_AFF4$p_value)

# NAやInfを含む行を削除
DF_cleaned_AFF4 <- DF_AFF4[!is.na(DF_AFF4$log2FC) & !is.infinite(DF_AFF4$log2FC) &
                             !is.na(DF_AFF4$p_value), ]
head(DF_cleaned_AFF4)
# 5. Volcano plotの作成
library(ggplot2)

# CBがTRUEの点を赤、FALSEの点を灰色で表示
volcano_AFF4 <- ggplot(DF_cleaned_AFF4, aes(x=log2FC, y=-log10(p_value))) +
  # 灰色の点を先に描画
  geom_point(data = subset(DF_cleaned_AFF4, CB == FALSE), aes(color = CB), size = 2) +
  # 赤い点を後に描画
  geom_point(data = subset(DF_cleaned_AFF4, CB == TRUE), aes(color = CB), size = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +  # 色の指定
  theme_minimal() +
  labs(x="log2 Fold Change", y="-log10(p-value)") +
  xlim(-5, 5) +   
  ylim(0, 7) +   
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14))
sum(DF_cleaned_AFF4$CB == TRUE)
# 6. 保存
ggsave("C:/Documents/20241117RNAseqIntergration/si7SK_AFF4_Volcano.pdf", plot = volcano_AFF4, width = 10, height = 8)
#################################################################################################################
# 必要なライブラリをロード
library(ggplot2)
# NAやInfを含む行を削除
DF_cleaned_AFF4 <- DF_AFF4[!is.na(DF_AFF4$log2FC) & !is.infinite(DF_AFF4$log2FC) &
                             !is.na(DF_AFF4$p_value), ]
head(DF_cleaned_AFF4)
# p-valueが0.05未満であるか、log2FCの絶対値が上位1000位以内のものをフィルタリング
DF_cleaned_AFF4$log_p_value <- -log10(DF_cleaned_AFF4$p_value)
DF_cleaned_AFF4$abs_log2FC <- abs(DF_cleaned_AFF4$log2FC)

# log2FCの絶対値で上位1000位を取得 (全行数が1000未満の場合はすべて含まれます)
top_1000_abs_log2FC <- sort(DF_cleaned_AFF4$abs_log2FC, decreasing = TRUE)[min(1000, nrow(DF_cleaned_AFF4))]

# 条件を満たすかどうかを判定
DF_cleaned_AFF4$highlight <- ifelse(DF_cleaned_AFF4$p_value < 0.05 & DF_cleaned_AFF4$abs_log2FC >= top_1000_abs_log2FC, "Significant", "Not Significant")

# 遺伝子名がRN7SKのデータポイントを特定してハイライト
DF_cleaned_AFF4$color <- ifelse(DF_cleaned_AFF4$external_gene_name == "RN7SK", "RN7SK", DF_cleaned_AFF4$highlight)

# RN7SKのデータポイントを抽出
RN7SK_row <- DF_cleaned_AFF4[DF_cleaned_AFF4$external_gene_name == "RN7SK", ]

# Volcano Plotを描画
# RN7SKのデータポイントを抽出
RN7SK_row <- DF_cleaned_AFF4[DF_cleaned_AFF4$external_gene_name == "RN7SK", ]
other_rows <- DF_cleaned_AFF4[DF_cleaned_AFF4$external_gene_name != "RN7SK", ]

# Volcano Plotを描画
volcano_plot_AFF4 <- ggplot() +
  # 他の点を描画 (最初に描画する)
  geom_point(data = other_rows, aes(x = log2FC, y = log_p_value, color = color), alpha = 1) +
  # RN7SKの点を描画 (最前面に描画)
  geom_point(data = RN7SK_row, aes(x = log2FC, y = log_p_value), color = "blue", size = 3, alpha = 1) +
  # RN7SKにラベルと線を追加
  geom_text(data = RN7SK_row, aes(x = log2FC, y = log_p_value, label = external_gene_name), 
            hjust = -0.2, vjust = 1.2, size = 5, color = "blue") +
  geom_segment(data = RN7SK_row, 
               aes(x = log2FC, y = log_p_value, xend = log2FC + 0.5, yend = log_p_value - 0.5), 
               arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
  scale_color_manual(values = c("Not Significant" = "black", "Significant" = "red")) +
  labs(
    title = "AFF4 RNAseq p_value < 0.05 & top_1000_abs_log2FC genes",
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  ) +
  xlim(-8, 8) +   
  ylim(0, 8) +   
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

# PDFに保存
output_path <- "C:/Documents/20241117RNAseqIntergration/volcano_plot_AFF4_KD_p_value005&top_1000_abs_log2FCgenes_with_label.pdf"
pdf(output_path, width = 8, height = 6) # サイズは適宜調整
print(volcano_plot_AFF4)
dev.off()

# 保存完了メッセージ
cat("PDFファイルが以下に保存されました:", output_path, "\n")
#################################################################################################################

# CB.1がTRUEの行を抽出
df_true <- DF_AFF4[DF_AFF4$CB.1 == TRUE, ]

# 行数をカウント
nrow(df_true)
# CB.1がTRUEで、p_valueが0.05以下、log2FCの絶対値が0.5より大きい行を抽出
df_filtered <- df_true[df_true$p_value <= 0.05 & abs(df_true$log2FC) > 0.5, ]
df_filtered_total_AFF4 <- DF_AFF4[DF_AFF4$p_value <= 0.05 & abs(DF_AFF4$log2FC) > 0.5, ]
# 抽出した行の数をカウント
nrow(df_filtered)
nrow(df_filtered_total_AFF4)

# 必要なパッケージを読み込む
library(ggplot2)

# Load necessary package
library(ggplot2)

# Create a data frame
data <- data.frame(
  category = c("Affected", "Not Affected"),
  count = c(56, 71)
)

# Calculate fractions
data$fraction = data$count / sum(data$count)

# Create labels for percentages
data$label = paste0(round(data$fraction * 100, 1), "%")

# Create pie chart
p <- ggplot(data, aes(x = "", y = count, fill = factor(category))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +  # Convert to pie chart
  theme_void() +  # Remove axes and background
  labs(title = "Proportion of Affected and Not Affected Genes") +
  scale_fill_manual(values = c("Affected" = "#66b3ff", "Not Affected" = "#99ff99")) +  # Color specification
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))  # Add percentage labels

# Save as PDF
ggsave("C:/Documents/20241117RNAseqIntergration/affected_gene_ratio_AFF4.pdf", plot = p, device = "pdf")

# Create a data frame
data_2 <- data.frame(
  category = c("CB_Affected", "Other_Affected"),
  count = c(56, 49390)
)

# Calculate fractions
data_2$fraction = data_2$count / sum(data_2$count)

# Create labels for percentages
data_2$label = paste0(round(data_2$fraction * 100, 1), "%")

# Create pie chart
p <- ggplot(data_2, aes(x = "", y = count, fill = factor(category))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +  # Convert to pie chart
  theme_void() +  # Remove axes and background
  labs(title = "AFF4 Proportion of Affected and Not Affected Genes") +
  scale_fill_manual(values = c("CB_Affected" = "#66b3ff", "Other_Affected" = "#99ff99")) +  # Color specification
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))  # Add percentage labels

# Save as PDF
ggsave("C:/Documents/20241117RNAseqIntergration/CB_Others_gene_ratio_AFF4.pdf", plot = p, device = "pdf")


###################################################################################################
# CB.1がTRUEの行を抽出
df_true_7SK <- DF_7SK[DF_7SK$CB.1 == TRUE, ]

# 行数をカウント
nrow(DF_7SK_CB)
# CB.1がTRUEで、p_valueが0.05以下、log2FCの絶対値が0.5より大きい行を抽出
df_filtered_7SK <- DF_7SK_CB[DF_7SK_CB$p_value <= 0.05 & abs(DF_7SK_CB$log2FC) > 0.5, ]
df_filtered_total_7SK <- DF_7SK[DF_7SK$p_value <= 0.05 & abs(DF_7SK$log2FC) > 0.5, ]

# 抽出した行の数をカウント
nrow(df_filtered_7SK)
nrow(df_filtered_total_7SK)
# Create a data frame
data_7SK_P <- data.frame(
  category = c("Affected", "Not Affected"),
  count = c(7, 55)
)

# Calculate fractions
data_7SK_P$fraction = data_7SK_P$count / sum(data_7SK_P$count)

# Create labels for percentages
data_7SK_P$label = paste0(round(data_7SK_P$fraction * 100, 1), "%")

# Create pie chart
p_7SK_P <- ggplot(data_7SK_P, aes(x = "", y = count, fill = factor(category))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +  # Convert to pie chart
  theme_void() +  # Remove axes and background
  labs(title = "Proportion of Affected and Not Affected Genes") +
  scale_fill_manual(values = c("Affected" = "#66b3ff", "Not Affected" = "#99ff99")) +  # Color specification
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))  # Add percentage labels

# Save as PDF
ggsave("C:/Documents/20241117RNAseqIntergration/affected_gene_ratio_7SK.pdf", plot = p_7SK_P, device = "pdf")
# Create a data frame
data_3 <- data.frame(
  category = c("CB_Affected", "Other_Affected"),
  count = c(7, 47869)
)

# Calculate fractionss
data_3$fraction = data_3$count / sum(data_2$count)

# Create labels for percentages
data_3$label = paste0(round(data_3$fraction * 100, 1), "%")

# Create pie chart
p <- ggplot(data_3, aes(x = "", y = count, fill = factor(category))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +  # Convert to pie chart
  theme_void() +  # Remove axes and background
  labs(title = "7SK Proportion of Affected and Not Affected Genes") +
  scale_fill_manual(values = c("CB_Affected" = "#66b3ff", "Other_Affected" = "#99ff99")) +  # Color specification
  geom_text(aes(label = label), position = position_stack(vjust = 0.5))  # Add percentage labels

# Save as PDF
ggsave("C:/Documents/20241117RNAseqIntergration/CB_Others_gene_ratio_7SK.pdf", plot = p, device = "pdf")

###################################################################################################
# CB = TRUEのデータのみを抽出
DF_AFF4_CB <- DF_AFF4[DF_AFF4$CB == TRUE, ]
# log2FCの計算
DF_AFF4_CB$log2FC <- log2(DF_AFF4_CB$AFF4_mean / DF_AFF4_CB$ICE1_mean)

# p値の計算 (t検定)
p_values <- apply(DF_AFF4_CB[, c("AFF4_mean", "ICE1_mean")], 1, function(x) {
  tryCatch({
    t.test(x[1], x[2])$p.value  # 例として2つのグループの平均を比較
  }, error = function(e) {
    return(NA)  # エラーが発生した場合はNAを返す
  })
})
DF_AFF4_CB$p_value <- p_values

# -log10(p-value)を計算して、Volcano plotに使用する
DF_AFF4_CB$log_p_value <- -log10(DF_AFF4_CB$p_value)

# NAやInfを含む行を削除
DF_AFF4_CB_cleaned <- DF_AFF4_CB[!is.na(DF_AFF4_CB$log2FC) & !is.infinite(DF_AFF4_CB$log2FC) &
                                   !is.na(DF_AFF4_CB$p_value), ]


# Volcano plotの作成
library(ggplot2)

ggplot(DF_AFF4_CB, aes(x = log2FC, y = -log10(p_value))) +
  geom_point(aes(color = p_value < 0.05 & abs(log2FC) > 0.5), size = 2) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2FC", y = "-log10(p-value)") +
  theme(legend.position = "none")

#################################################################################################
