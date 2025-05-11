library(umap)
library(ggplot2)
library(dplyr)

# データフレームの読み込み（仮定）
df <- read.csv("Mean_summary_with_averages.csv")
setwd("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed")
# 分散計算: 各サンプルの変動（std）を計算
variance <- apply(df[, c("input_Ave", "Pol2_Ave", "Ser2_Ave", "Ser5_Ave")], 1, var)

# 分散上位10％の遺伝子を抽出
threshold <- quantile(variance, 0.9)
top_genes <- df[variance >= threshold, ]

# UMAPの適用
umap_model <- umap(top_genes[, c("input_Ave", "Pol2_Ave", "Ser2_Ave", "Ser5_Ave")])

# UMAP結果をデータフレーム化
umap_data <- data.frame(umap_model$layout)

# 遺伝子タイプをグループ化: 特定のRNAタイプを"others"にまとめる
valid_types <- c("lincRNA","protein_coding", "miRNA", "misc_RNA", "ribozyme", "rRNA", "scaRNA", "snoRNA", "snRNA", "sRNA", "TEC", "vault_RNA")

top_genes$group <- ifelse(
  top_genes$gene_biotype %in% valid_types,  # 上記のRNAタイプはそのまま表示
  top_genes$gene_biotype,  # 一致するものはそのまま
  "others"  # 一致しないものは"others"
)

# UMAPの結果に遺伝子タイプの情報を追加
umap_data$group <- top_genes$group

# UMAPの結果をデータフレーム化
umap_data <- data.frame(umap_model$layout)

# UMAPの結果の列名を確認（V1, V2がない場合に備えて）
colnames(umap_data)
colnames(umap_data) <- c("UMAP1", "UMAP2")# UMAPの結果がV1, V2でない場合、列名を変更する
umap_data$group <- top_genes$group# 選択した遺伝子のUMAP座標を確認

# UMAPプロット
ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point(alpha = 0.5, size = 1.5) +  # 点を小さめにして視認性を確保
  scale_color_manual(values = c("others" = "grey", 
                                "lincRNA" = "blue", 
                                "miRNA" = "green", 
                                "misc_RNA" = "purple", 
                                "ribozyme" = "orange", 
                                "rRNA" = "pink", 
                                "scaRNA" = "brown", 
                                "snoRNA" = "yellow", 
                                "snRNA" = "red", 
                                "sRNA" = "magenta", 
                                "TEC" = "black",
                                "protein_coding"="cyan",
                                "vault_RNA" = "darkgreen")) +
  theme_minimal() +
  labs(title = "UMAP of top 10% variance genes", x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "right")

# UMAPの条件に一致するデータを抽出
selected_umap_data <- umap_data[umap_data$UMAP1 < -10 & umap_data$UMAP2 > 12, ]
head(selected_umap_data)
head(df)
# selected_umap_dataのインデックスを使って、遺伝子名をdfから取得
selected_gene_ids <- rownames(selected_umap_data)

# dfから対応する遺伝子名を取得
selected_genes_info <- df[selected_gene_ids, c("ENSG","external_gene_name", "gene_biotype")]

# 結果を表示
selected_genes_info

snrna_count <- sum(selected_genes_info$gene_biotype == "snRNA")# snRNAタイプの遺伝子数をカウント
total_count <- nrow(selected_genes_info)# 全遺伝子数をカウント
snrna_percentage <- (snrna_count / total_count) * 100# snRNAの割合を計算
snrna_percentage# 結果を表示

# selected_genes_info を CSV ファイルとして保存
write.csv(selected_genes_info, file = "C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/selected_genes_info.csv", row.names = FALSE)
#########################################################################################################################################################
# 必要なライブラリの読み込み
library(ggplot2)
library(Rtsne)
library(umap)

# 分散計算: 各サンプルの変動（std）を計算
variance <- apply(df[, c("input_Ave", "Pol2_Ave", "Ser2_Ave", "Ser5_Ave")], 1, var)

# 分散上位10％の遺伝子を抽出
threshold <- quantile(variance, 0.9)
top_genes <- df[variance >= threshold, ]

# ここではUMAPを使った結果を再確認
umap_model <- umap(top_genes[, c("input_Ave", "Pol2_Ave", "Ser2_Ave", "Ser5_Ave")])
umap_data <- data.frame(umap_model$layout)
colnames(umap_data) <- c("UMAP1", "UMAP2")

# t-SNEの適用
tsne_model <- Rtsne(top_genes[, c("input_Ave", "Pol2_Ave", "Ser2_Ave", "Ser5_Ave")])
tsne_data <- data.frame(tsne_model$Y)
colnames(tsne_data) <- c("tSNE1", "tSNE2")

# 重複行を削除
top_genes_no_duplicates <- top_genes[!duplicated(top_genes[, c("input_Ave", "Pol2_Ave", "Ser2_Ave", "Ser5_Ave")]), ]

# t-SNEを再実行
tsne_model <- Rtsne(top_genes_no_duplicates[, c("input_Ave", "Pol2_Ave", "Ser2_Ave", "Ser5_Ave")])
tsne_data <- data.frame(tsne_model$Y)
colnames(tsne_data) <- c("tSNE1", "tSNE2")

# PCAの適用
pca_model <- prcomp(top_genes[, c("input_Ave", "Pol2_Ave", "Ser2_Ave", "Ser5_Ave")], scale = TRUE)
pca_data <- data.frame(pca_model$x[, 1:2])  # PC1とPC2を取得
colnames(pca_data) <- c("PCA1", "PCA2")

# 遺伝子タイプをグループ化: 特定のRNAタイプを"others"にまとめる
valid_types <- c("lincRNA", "miRNA", "misc_RNA", "ribozyme", "rRNA", "scaRNA", "snoRNA", "snRNA", "sRNA", "TEC", "vault_RNA")

top_genes$group <- ifelse(
  top_genes$gene_biotype %in% valid_types,  # 上記のRNAタイプはそのまま表示
  top_genes$gene_biotype,  # 一致するものはそのまま
  "others"  # 一致しないものは"others"
)

# グループ情報を各プロットに追加
umap_data$group <- top_genes$group
# 重複を削除したデータのグループ情報を保持
tsne_data$group <- top_genes_no_duplicates$group
pca_data$group <- top_genes$group
# group列が追加されたか確認
tsne_data$group <- top_genes_no_duplicates$group
# top_genes_no_duplicatesにgroup列があるか確認
head(top_genes_no_duplicates)

head(tsne_data)


library(ggplot2)
# ggplot2のバージョン確認
packageVersion("ggplot2")

# UMAPプロット
umap_plot <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("others" = "grey", 
                                "lincRNA" = "blue", 
                                "miRNA" = "green", 
                                "misc_RNA" = "purple", 
                                "ribozyme" = "orange", 
                                "rRNA" = "pink", 
                                "scaRNA" = "brown", 
                                "snoRNA" = "yellow", 
                                "snRNA" = "red", 
                                "sRNA" = "magenta", 
                                "TEC" = "grey", 
                                "vault_RNA" = "darkgreen")) +
  theme_minimal() +
  labs(title = "UMAP of top 10% variance genes", x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "right")
# UMAPプロットの表示
print(umap_plot)

# t-SNEプロット

tsne_plot <- ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, color = group)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("others" = "grey", 
                                "lincRNA" = "blue", 
                                "miRNA" = "green", 
                                "misc_RNA" = "purple", 
                                "ribozyme" = "orange", 
                                "rRNA" = "pink", 
                                "scaRNA" = "brown", 
                                "snoRNA" = "yellow", 
                                "snRNA" = "red", 
                                "sRNA" = "magenta", 
                                "TEC" = "grey", 
                                "vault_RNA" = "darkgreen")) +
  theme_minimal() +
  labs(title = "t-SNE of top 10% variance genes", x = "tSNE1", y = "tSNE2") +
  theme(legend.position = "right")

# プロットを表示
print(tsne_plot)
# PCAプロット
pca_plot <- ggplot(pca_data, aes(x = PCA1, y = PCA2, color = group)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("others" = "grey", 
                                "lincRNA" = "blue", 
                                "miRNA" = "green", 
                                "misc_RNA" = "purple", 
                                "ribozyme" = "orange", 
                                "rRNA" = "pink", 
                                "scaRNA" = "brown", 
                                "snoRNA" = "yellow", 
                                "snRNA" = "red", 
                                "sRNA" = "magenta", 
                                "TEC" = "grey", 
                                "vault_RNA" = "darkgreen")) +
  theme_minimal() +
  labs(title = "PCA of top 10% variance genes", x = "PCA1", y = "PCA2") +
  theme(legend.position = "right")


install.packages("gridExtra")
# プロットの表示
library(gridExtra)
grid.arrange(umap_plot, tsne_plot, pca_plot, ncol = 3)

# 特定のRNAタイプを"valid_types"として設定
valid_types <- c("lincRNA", "miRNA", "misc_RNA", "ribozyme", "rRNA", 
                 "scaRNA", "snoRNA", "snRNA", "sRNA", "TEC", "vault_RNA")

# gene_biotypeがvalid_typesに一致する場合はそのまま、そうでない場合は"others"にする
top_genes_no_duplicates$group <- ifelse(
  top_genes_no_duplicates$gene_biotype %in% valid_types, 
  top_genes_no_duplicates$gene_biotype, 
  "others"
)

# group列が追加されたか確認
head(top_genes_no_duplicates)

