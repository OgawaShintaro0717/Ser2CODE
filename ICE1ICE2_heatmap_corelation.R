"C:\Documents\ICE2_ICE1_Pol2_Coilin\average_signal"
# R script to calculate Spearman correlations and create heatmap
library(pheatmap)
library(readr)
install.packages("readr")

# ディレクトリ内のファイル指定
signal_files <- list.files(path = "C:/Documents/ICE2_ICE1_Pol2_Coilin/average_signal/", pattern = "*.signal.txt", full.names = TRUE)

# 各ファイルをデータフレームに読み込んで、統一してマージします
signal_data <- lapply(signal_files, function(file) {
  sample <- tools::file_path_sans_ext(basename(file))  # ファイル名からサンプル名取得
  df <- read_delim(file, delim = "\t")
  df$sample <- sample
  return(df)
})

df_combined <- do.call(rbind, signal_data)

# 各データフレームの列名を表示する
lapply(signal_data, colnames)
# 列名を統一する関数
signal_data <- lapply(signal_data, function(df) {
  # 列名を適切に設定する
  colnames(df) <- c("gene_id", "expression1", "expression2", "value1", "value2", "value3", "sample")
  return(df)
})

# 結合する
df_combined <- do.call(rbind, signal_data)
head(df_combined)
# 結果を出力する
write.csv(df_combined, 'merged_signal_data.csv', row.names = FALSE)

cat('Signal data merged and saved successfully!')

# パッケージのインストール（必要に応じて）
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("corrplot")) install.packages("corrplot")

# ライブラリの読み込み
library(pheatmap)

# サンプルごとの遺伝子発現データだけを抽出
df_corr <- df_combined[, -c(1, 7)]  # gene_idとsample列以外を選択

# データフレームを行ごとの遺伝子、列ごとのサンプルに整形する
df_matrix <- as.matrix(df_corr)

# スピアマン相関係数の計算
cor_matrix <- cor(df_matrix, method = "spearman")

# ヒートマップの作成
pheatmap(cor_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Spearman Correlation Heatmap",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_legend = TRUE
)


