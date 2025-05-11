install.packages("DESeq2") 
library("DESeq2")
library(ggplot2)

# CSVファイルの読み込み
data <- read.csv("C:/Documents/241117ICE_AFF4_RNAseq/expected_count_matome.csv", row.names = 1)
# データを確認
head(data)
# サンプルの条件を作成
coldata <- data.frame(
  condition = factor(rep(c("group1", "group2"), each = 3)),
  row.names = colnames(data)
)
data_int <- round(data)

# DESeq2のデータセットを作成
dds <- DESeqDataSetFromMatrix(countData = data_int, colData = coldata, design = ~ condition)

# DESeq2解析を実行
dds <- DESeq(dds)

# 結果を取得
res <- results(dds)
head(res)
res$negLogPValue <- -log10(res$pvalue) # -log10(p-value)を計算

# 上限を設定し、それを超える点を別の形や色で表示
res$negLogPValueCapped <- ifelse(res$negLogPValue > 300, 300, res$negLogPValue)  # y軸の最大値でキャップ
res$isCapped <- res$negLogPValue > 300  # 上限を超えたかどうかをフラグにする

ggplot(res, aes(x = log2FoldChange, y = negLogPValueCapped)) +
  geom_point(aes(color = pvalue < 0.05 & abs(log2FoldChange) > 1), size = 2) +
  geom_point(data = subset(res, isCapped), shape = 24, size = 3, fill = "black") +  # キャップした点を三角形で表示
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  xlim(c(-10, 10)) +
  ylim(c(0, 300)) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(p-value)") +
  theme(legend.position = "none")



# 差次的発現した遺伝子をフィルタリング
significant_genes <- subset(res, pvalue < 0.05 & log2FoldChange > 1)
head(significant_genes)
# 遺伝子名（rownames）を取得
gene_list <- rownames(significant_genes)

# 結果を確認
print(gene_list)
# biomaRtパッケージのインストールとロード
library(biomaRt)

# Ensemblのデータベースに接続
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 遺伝子IDから遺伝子シンボルを取得
gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(significant_genes),
  mart = mart
)

# 遺伝子IDを列として追加
significant_genes <- data.frame(
  ensembl_gene_id = rownames(significant_genes),
  significant_genes
)

# 結果をデータフレームにマージ
significant_genes <- data.frame(
  ensembl_gene_id = rownames(significant_genes),
  significant_genes
)
significant_genes <- merge(significant_genes, gene_symbols, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)

# 結果を確認
head(significant_genes)

write.csv(significant_genes, "C:/Documents/241117ICE_AFF4_RNAseq/Increasing_significant_genes_with_symbols.csv", row.names = FALSE)


# 減少した遺伝子をフィルタリング
downregulated_genes <- subset(res, pvalue < 0.05 & log2FoldChange < -1)

# 遺伝子IDを列として追加
downregulated_genes <- data.frame(
  ensembl_gene_id = rownames(downregulated_genes),
  downregulated_genes
)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(res),
  mart = mart
)

# 一般名をマージ
downregulated_genes <- merge(
  downregulated_genes,
  gene_symbols,
  by.x = "ensembl_gene_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)
head(downregulated_genes)
 
# 結果をファイルに保存
write.csv(downregulated_genes, "C:/Documents/241117ICE_AFF4_RNAseq/downregulated_genes_with_symbols.csv", row.names = FALSE)



# PCAを実行して結果を可視化
vsd <- vst(dds, blind = FALSE)  # VST (Variance Stabilizing Transformation)を適用
pca_result <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)

# PCA結果を表示
head(pca_result)

# PCAの実行
vsd <- vst(dds, blind = FALSE)  # VST変換
pca_res <- prcomp(t(assay(vsd)))  # 転置してPCAを実行
percent_var <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100 # 分散の割合を計算

# PCA結果をデータフレームにまとめる
pca_result <- data.frame(pca_res$x)
pca_result$condition <- coldata$condition  # condition列を追加（coldataが適切であれば）
pca_result$sample <- rownames(pca_result)  # サンプル名を追加

# 分散割合を追加
pca_result$PC1_percentVar <- percent_var[1]
pca_result$PC2_percentVar <- percent_var[2]
# PCAの可視化
ggplot(pca_result, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste("PC1: ", round(pca_result$PC1_percentVar, 1), "% variance")) +
  ylab(paste("PC2: ", round(pca_result$PC2_percentVar, 1), "% variance")) +
  ggtitle("PCA of RNA-Seq data") +
  theme_minimal() +
  theme(legend.position = "top")

ggplot(data, aes(x = condition, y = expression)) +
  geom_boxplot() +
  theme_minimal()

head(data)
data$condition <- rep(c("group1", "group2"), each = 4)  # 必要に応じて修正

# パッケージのインストール
install.packages("Rtsne")
library(Rtsne)

# t-SNEの実行（初期データは通常正規化データを推奨）
tsne_result <- Rtsne(t(data), dims = 2, perplexity = 1)

# 結果をデータフレーム化
tsne_df <- as.data.frame(tsne_result$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")
tsne_df$condition <- coldata$condition

# プロット
ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "t-SNE Plot", x = "t-SNE 1", y = "t-SNE 2")
print(coldata)

# パッケージのインストール
install.packages("umap")
library(umap)

# デフォルト設定を変更
custom_config <- umap.defaults
custom_config$n_neighbors <- 2  # サンプル数（6）未満の値に設定

# UMAP 実行
umap_result <- umap(t(data), config = custom_config)

# 結果をデータフレーム化
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$condition <- coldata$condition

# プロット
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2")



library(tidyr)
# tibble パッケージのインストール
install.packages("tibble")

# tibble パッケージのロード
library(tibble)
# 長い形式に変換
long_data <- data %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(cols = starts_with("sample"), names_to = "sample", values_to = "expression")

# Boxplot の作成
ggplot(long_data, aes(x = sample, y = expression)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Gene Expression Across Samples",
       x = "Sample",
       y = "Expression Level")
library(dplyr)

# 四分位範囲の計算
iqr_values <- long_data %>%
  summarize(
    lower = quantile(expression, 0.25),
    upper = quantile(expression, 0.75),
    iqr = IQR(expression)
  )

# 上限の計算（75%分位点 + IQRの1.5倍）
y_max <- iqr_values$upper + 1.5 * iqr_values$iqr

# Boxplot の作成
ggplot(long_data, aes(x = sample, y = expression)) +
  geom_boxplot() +
  theme_minimal() +
  ylim(0, y_max) +  # y軸の上限を設定
  labs(title = "Gene Expression Across Samples",
       x = "Sample",
       y = "Expression Level")

ggplot(long_data, aes(x = sample, y = expression)) +
  geom_boxplot(coef = 1) +  # IQRの1倍に変更
  theme_minimal() +
  labs(title = "Gene Expression Across Samples",
       x = "Sample",
       y = "Expression Level")


install.packages("pheatmap")
library(pheatmap)
pheatmap(assay(vsd), cluster_rows = TRUE, cluster_cols = TRUE, scale = "row")

# パッケージのインストール
install.packages("fastICA")
library(fastICA)

# ICAの実行
ica_result <- fastICA(t(data), n.comp = 2)

# 結果をデータフレーム化
ica_df <- as.data.frame(ica_result$S)
colnames(ica_df) <- c("IC1", "IC2")
ica_df$condition <- coldata$condition

# プロット
ggplot(ica_df, aes(x = IC1, y = IC2, color = condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "ICA Plot", x = "IC 1", y = "IC 2")

