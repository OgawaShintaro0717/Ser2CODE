ASCL1_6h_ChIP_1_1000output.csv
setwd("C:/Documents/all_1000bpcsv_files")
list.files("C:/Documents/all_1000bpcsv_files")
directory_path <- "C:/Documents/all_1000bpcsv_files"
filtered_files <- list.files(directory_path, pattern = "_1000output.csv$")
files_to_delete <- setdiff(list.files(directory_path), filtered_files)
for (file in files_to_delete) {
  file_path <- file.path(directory_path, file)
  file.remove(file_path)  # ファイルを削除
  cat("Deleted:", file_path, "\n")  # 削除したファイルのパスを表示
}
list.files("C:/Documents/all_1000bpcsv_files")
ASCL1_6h <- read.csv("ASCL1_6h_ChIP_1_1000output.csv")
head(ASCL1_6h)

# 対象ディレクトリのパス
directory_path <- "C:/Documents/all_1000bpcsv_files"

# すべての対象ファイルをリスト（_1000output.csv のファイル）
filtered_files <- list.files(directory_path, pattern = "_1000output.csv$")

# 最初のCSVファイルを読み込んで、ENSG列を取り出す（これを基準に他の列を追加）
first_file_path <- file.path(directory_path, filtered_files[1])
first_df <- read.csv(first_file_path)

# 空のデータフレームを作成、最初のファイルのENSG列を最初の列として追加
combined_data <- data.frame(ENSG = first_df$ENSG, stringsAsFactors = FALSE)

# 各ファイルの処理
for (file in filtered_files) {
  # 各CSVファイルを読み込む
  file_path <- file.path(directory_path, file)
  df <- read.csv(file_path)
  
  # ファイル名（拡張子除く）を列名として設定
  file_name <- sub("_1000output.csv$", "", file)
  
  # "Mean"列を抽出し、ファイル名を列名としてデータフレームに追加
  combined_data[[file_name]] <- df$Mean
}

# 結果をCSVファイルに保存
write.csv(combined_data, "combined_means.csv", row.names = FALSE)

# 結果を確認
head(combined_data)
必要なパッケージを読み込み
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(pheatmap)

combined_data <- read.csv("C:/Documents/all_1000bpcsv_files/combined_means.csv")
# データフレームから必要な列を選択
columns_of_interest <- c("N_81_221_6h_2_CPM_normalized", 
                         "N_81_272_6h_1_CPM_normalized", 
                         "N_81_272_6h_2_CPM_normalized", 
                         "N_81_190_6h_1_CPM_normalized", 
                         "N_81_190_6h_2_CPM_normalized", 
                         "N_81_221_6h_1_CPM_normalized", 
                         "N_81_179_6h_1_CPM_normalized", 
                         "N_81_179_6h_2_CPM_normalized", 
                         "GFP_6h_ChIP_1", 
                         "GFP_6h_ChIP_2")

columns_of_interest <- c("ASCL1_6h_ChIP_1", 
                         "ASCL1_6h_ChIP_2", 
                         "MYOD1_6h_ChIP_1", 
                         "MYOD1_6h_ChIP_2", 
                         "NGN2_6h_ChIP_1", 
                         "NGN2_6h_ChIP_2", 
                         "GFP_6h_ChIP_1", 
                         "GFP_6h_ChIP_2")

# 必要な列を抽出
selected_data <- combined_data[, columns_of_interest]
head(selected_data)
N_81_221_r2 N_81_272_r1 N_81_272_r2 N_81_190_r1 N_81_190_r2 N_81_221_r1 N_81_179_r1 N_81_179_r2 GFP_r1 GFP_r
new_column_names <- c("N_81_221_r2", 
                      "N_81_272_r1", 
                      "N_81_272_r2", 
                      "N_81_190_r1", 
                      "N_81_190_r2", 
                      "N_81_221_r1", 
                      "N_81_179_r1", 
                      "N_81_179_r2", 
                      "GFP_r1", 
                      "GFP_r2")
colnames(selected_data) <- new_column_names
head(selected_data)
# スピアマン相関行列を計算
cor_matrix <- cor(selected_data, method = "spearman")

# 相関行列をヒートマップで表示
cor_matrix_melted <- melt(cor_matrix)  # 相関行列を長い形式に変換
# カスタムカラーパレット
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# ブレークポイントを調整
custom_breaks <- c(seq(-1, 0.2, length.out = 50), seq(0.21, 1, length.out = 50))
pheatmap(cor_matrix,
         color = custom_colors,
         breaks = custom_breaks,
         main = "Spearman Correlation Heatmap",
         fontsize_row = 10,
         fontsize_col = 10,
         cellwidth = 10,       # セルの幅を指定（適宜調整）
         cellheight = 10,
         display_numbers = F# セルの高さを指定（幅と同じ値にする）
)
#############################################################################
# "N_81_179_r2" 列を除外
selected_data <- selected_data[, !colnames(selected_data) %in% "N_81_179_r2"]

# スピアマン相関行列を再度計算
cor_matrix <- cor(selected_data, method = "spearman")

# カスタムカラーパレット
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# ブレークポイントを調整
custom_breaks <- c(seq(-1, 0.2, length.out = 50), seq(0.21, 1, length.out = 50))
pheatmap(cor_matrix,
         color = custom_colors,
         breaks = custom_breaks,
         main = "Spearman Correlation Heatmap",
         fontsize_row = 10,
         fontsize_col = 10,
         cellwidth = 10,       # セルの幅を指定（適宜調整）
         cellheight = 10,
         display_numbers = F)  # セルの高さを指定（幅と同じ値にする）
# "N_81_179_r2", "N_81_190_r1", "N_81_190_r2" 列を除外
selected_data <- selected_data[, !colnames(selected_data) %in% c("N_81_179_r2", "N_81_190_r1", "N_81_190_r2")]

# スピアマン相関行列を再度計算
cor_matrix <- cor(selected_data, method = "spearman")

# カスタムカラーパレット
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# ブレークポイントを調整
custom_breaks <- c(seq(-1, 0.2, length.out = 50), seq(0.21, 1, length.out = 50))
pheatmap(cor_matrix,
         color = custom_colors,
         breaks = custom_breaks,
         main = "Spearman Correlation Heatmap ",
         fontsize_row = 10,
         fontsize_col = 10,
         cellwidth = 10,       # セルの幅を指定（適宜調整）
         cellheight = 10,
         display_numbers = F)  # セルの高さを指定（幅と同じ値にする）

#############################################################################
head(selected_data)
# データの転置
selected_data_t <- t(selected_data_no_zero)
selected_data_no_zero_t <- selected_data_t[, colSums(selected_data_t != 0) > 0]

# PCAの実行
pca <- prcomp(selected_data_no_zero_t, scale. = TRUE)

# PCA結果をデータフレームに変換
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])

# PCAプロット
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = rownames(pca_df))) +  # サンプル名で色付け
  theme_minimal() +
  labs(title = "PCA of Samples (Excluding Zero-Only Columns)", x = "Principal Component 1", y = "Principal Component 2")

####################################################################
library(umap)
library(ggplot2)

# 2. UMAPの実行
umap_model <- umap(selected_data_no_zero_t, n_neighbors = min(nrow(selected_data_no_zero_t) - 1, 15))

# 3. UMAP結果をデータフレームに変換
umap_df <- data.frame(UMAP1 = umap_model$layout[, 1], UMAP2 = umap_model$layout[, 2])

# 4. UMAPプロット
ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = rownames(umap_df)), size = 4) +  # サンプル名で色付け
  theme_minimal() +
  labs(title = "UMAP of Samples (Excluding Zero-Only Columns)", x = "UMAP1", y = "UMAP2")
#############################################################################
# k-meansクラスタリング（クラスタ数を決める）
set.seed(42)  # 再現性のために乱数の種を設定
kmeans_result <- kmeans(selected_data_no_zero_t, centers = 2)  # 例えば3クラスタに分ける

# UMAPの結果にクラスタラベルを追加
umap_df$cluster <- as.factor(kmeans_result$cluster)

# UMAPプロット
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 4) +  # 点の大きさを変更
  theme_minimal() +
  labs(title = "UMAP of Samples with k-means Clustering", x = "UMAP1", y = "UMAP2")
############################################################################################
library(Rtsne)
tsne_results <- Rtsne(selected_data_no_zero_t, dims = 2, perplexity = 3)
tsne_df <- data.frame(TSNE1 = tsne_results$Y[, 1], TSNE2 = tsne_results$Y[, 2], Sample = rownames(selected_data_no_zero_t))

ggplot(tsne_df, aes(x = TSNE1, y = TSNE2)) +
  geom_point(aes(color = Sample), size = 4) +
  theme_minimal() +
  labs(title = "t-SNE of Samples", x = "t-SNE1", y = "t-SNE2")
###########################################################################################
distance_matrix <- dist(selected_data_no_zero_t, method = "euclidean")  # 距離行列
hc <- hclust(distance_matrix, method = "ward.D2")  # 階層型クラスタリング

# デンドログラムのプロット
plot(hc, labels = rownames(selected_data_no_zero_t), main = "Hierarchical Clustering")
############################################################################################
##Cosine Similarity
install.packages("lsa")
library(lsa)

# 指定された新しい列名
new_colnames <- c(
  "ASCL1_r1", "ASCL1_r2", "MYOD1_r1", "MYOD1_r2", 
  "NGN2_r1", "NGN2_r2","GFP_r1", "GFP_r2"
)

rownames(selected_data) <- new_colnames

# コサイン類似度の計算（既存のコードを利用）
cosine_similarity <- cosine(t(selected_data_no_zero_t))
head(selected_data)

# ヒートマップのプロット
pheatmap(
  cosine_similarity,
  main = "Cosine Similarity Heatmap",
  labels_col = new_colnames,  # 列名をラベルとして設定
  labels_row = new_colnames,  # 行名をラベルとして設定
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,  # セルの幅を指定（適宜調整）
  cellheight = 15  # セルの高さを指定（幅と一致させる）
)

# N_81_179_r2 を除外
selected_data_filtered <- selected_data_no_zero_t[new_colnames != "N_81_179_r2", ]
filtered_colnames <- new_colnames[new_colnames != "N_81_179_r2"]

# コサイン類似度の計算
cosine_similarity <- cosine(t(selected_data_filtered))

# ヒートマップのプロット
pheatmap(
  cosine_similarity,
  main = "Cosine Similarity Heatmap (Excluding N_81_179_r2)",
  labels_col = filtered_colnames,  # 除外後の列名をラベルとして設定
  labels_row = filtered_colnames,  # 除外後の行名をラベルとして設定
  fontsize_row = 10,
  fontsize_col = 10,
  cellwidth = 15,  # セルの幅を指定（適宜調整）
  cellheight = 15  # セルの高さを指定（幅と一致させる）
)
##################################################################################################
mds <- cmdscale(distance_matrix, k = 2)
mds_df <- data.frame(MDS1 = mds[, 1], MDS2 = mds[, 2], Sample = rownames(selected_data_no_zero_t))

ggplot(mds_df, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = Sample), size = 4) +
  theme_minimal() +
  labs(title = "MDS of Samples", x = "MDS1", y = "MDS2")
############################################################################################

# ヒートマップのプロット
pheatmap(cosine_similarity, main = "Cosine Similarity Heatmap")

# データの読み込み
data <- read.csv("C:/Documents/all_csv_files/combined_mean_columns.csv")

# ENSG列以外のデータ（サンプルごとのリードカウント行列）
expression_data <- data[, -1]

# Spearman相関行列の計算
cor_matrix <- cor(expression_data, method = "spearman")

# 相関プロットの作成
library(pheatmap)

# 相関行列をヒートマップで可視化
pheatmap(cor_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Spearman Correlation Heatmap",
         fontsize_row = 10,   # 行名の文字サイズを設定
         fontsize_col = 10,   # 列名の文字サイズを設定
         show_rownames = FALSE,
         show_colnames = TRUE)

packageVersion("stats")

setwd("C:/Documents/all_csv_files")
# ディレクトリ内の全CSVファイルリストを取得
file_list <- list.files(pattern = "_output.csv$")
# 空のデータフレームを作成しておく
final_data <- data.frame(ENSG = character())

# 各CSVファイルをループして処理する
for (file in file_list) {
  # CSVファイルを読み込む
  data <- read.csv(file)
  
  # ファイル名から列名を生成（例： ASCL1_6h_ChIP_1）
  sample_name <- gsub("_output.csv", "", file)
  
  # `ENSG`と`Mean`列のみ抽出
  mean_data <- data[, c("ENSG", "Mean")]
  colnames(mean_data) <- c("ENSG", sample_name)
  
  # 最終データフレームに統一して結合
  if (nrow(final_data) == 0) {
    final_data <- mean_data
  } else {
    final_data <- merge(final_data, mean_data, by = "ENSG")
  }
}

# 結果をCSVファイルに出力
write.csv(final_data, "combined_mean_data.csv", row.names = FALSE)
final_data <- read.csv("combined_mean_data.csv")
cat("データをcombined_mean_data.csvにまとめました！\n")
head(final_data)

# 必要なライブラリの読み込み
library(pheatmap)

# データの読み込み（すでにある`final_data`を使用）
# Spearman相関行列の計算
expression_data <- final_data[, -1]  # ENSG列を除く
cor_matrix <- cor(expression_data, method = "spearman")

# Spearman相関プロット作成（ヒートマップ）
pheatmap(cor_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Spearman Correlation Heatmap",
         fontsize_row = 10,
         fontsize_col = 10,
         show_rownames = TRUE,
         show_colnames = TRUE,
         display_numbers = T)
colnames(final_data)
selected_columns <- final_data[, c("N_81_221_6h_2_CPM_normalized",
                                   "N_81_190_6h_2_CPM_normalized",
                                   "N_81_179_6h_1_CPM_normalized",
                                   "N_81_272_6h_2_CPM_normalized",
                                   "MYOD1_6h_ChIP_2",
                                   "NGN2_6h_ChIP_1",
                                   "ASCL1_6h_ChIP_1", 
                                   "GFP_6h_ChIP_1")]
head(selected_columns)
colnames(selected_columns)
# 列名を変更する
colnames(selected_columns) <- c("N_81_221", "N_81_190", "N_81_179", "N_81_272","MYOD1","NGN2","ASCL1","GFP")

# Spearman相関行列の計算
expression_data <- selected_columns  # ENSG列を除く
head(expression_data)
cor_matrix <- cor(expression_data, method = "spearman")

# Spearman相関プロット作成（ヒートマップ）
pheatmap(cor_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Spearman Correlation Heatmap",
         fontsize_row = 5,
         fontsize_col = 5,
         show_rownames = TRUE,
         show_colnames = TRUE)
# カスタムカラーパレット
custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# ブレークポイントを調整
custom_breaks <- c(seq(0, 0.48, length.out = 50), seq(0.52, 1, length.out = 50))
pheatmap(cor_matrix,
         color = custom_colors,
         breaks = custom_breaks,
         main = "Spearman Correlation Heatmap",
         fontsize_row = 10,
         fontsize_col = 10,
         cellwidth = 20,       # セルの幅を指定（適宜調整）
         cellheight = 20,
         display_numbers = F# セルの高さを指定（幅と同じ値にする）
)

#############################################################################
library(umap)
library(ggplot2)
# データフレームを転置
transposed_data <- t(selected_columns)
# UMAPで次元削減を実施
umap_result <- umap(as.matrix(selected_columns))

# 結果をデータフレームに変換
umap_df <- data.frame(UMAP1 = umap_result$layout[,1],
                      UMAP2 = umap_result$layout[,2])

# ggplotでプロット作成
ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(color = "steelblue", size = 3) +
  labs(title = "UMAP Plot of Selected Samples") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14))






