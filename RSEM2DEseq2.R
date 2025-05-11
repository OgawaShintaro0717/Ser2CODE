if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("biomaRt")
a

setwd("C:/Documents/241114si7SKRNAseq")#/SRR5818141.csv

# 必要なライブラリを読み込む
library(biomaRt)

# RSEMカウントデータの読み込み
# カウントデータには "Ensembl_ID" という列が含まれていると仮定

rsem_data <- read.csv("C:/Documents/241114si7SKRNAseq/SRR5818148.csv")  # ファイル名を適宜変更

# Ensemblのデータベースに接続
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # ヒトデータを例にしています

# Ensembl ID、Gene Symbol、Transcript Typeを取得
gene_info <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol", "transcript_biotype"),
  values = rsem_data$gene_id,
  mart = ensembl
)

# RSEMデータにGene SymbolとTranscript Typeをマージ
rsem_annotated <- merge(rsem_data, gene_info, by.x = "gene_id", by.y = "ensembl_gene_id", all.x = TRUE)

# 結果を保存
write.csv(rsem_annotated, "SRR5818148_rsem_counts_with_gene_symbols_and_transcript_type.csv", row.names = FALSE)


# データの読み込み
data <- read.csv("C:/Documents/241114si7SKRNAseq/SRR5818148_rsem_counts_with_gene_symbols_and_transcript_type.csv")  # ファイル名を適切に変更

# 必要な列の抽出
selected_data <- data[, c("TPM", "hgnc_symbol")]

# 重複行の削除
unique_data <- unique(selected_data)

# 結果を確認
head(unique_data)

# CSVファイルに保存
write.csv(unique_data, file = "SRR5818148_filtered_data.csv", row.names = FALSE)


library("dplyr")
# ファイル名をリスト化（パターンを適切に変更）
library(dplyr)
library(purrr)

# ファイル名をリスト化
file_names <- list.files(pattern = "SRR581814[1-8]_filtered_data\\.csv")  # ファイル名パターンを指定

# 各ファイルをリストとして読み込む
data_list <- file_names %>%
  set_names() %>%
  map(~ read.csv(.))
# 最初のデータフレームを確認
head(data_list[[1]])
# データフレームを統合
merged_data <- data_list %>%
  reduce(function(df1, df2) {
    # 'hgnc_symbol'で結合し、TPM列の名前にサフィックスを追加
    full_join(df1, df2, by = "hgnc_symbol", suffix = c("", paste0("_", tools::file_path_sans_ext(basename(names(df2)))))) 
  })

merged_data <- data_list %>%
  reduce(function(df1, df2) {
    # ファイル名を適切に取得し、サフィックスを2つずつ作成
    file_suffix <- tools::file_path_sans_ext(basename(names(df2)))
    full_join(df1, df2, by = "hgnc_symbol", suffix = c("", file_suffix))
  })






# ファイルのパスを指定
file_path <- "C:/Documents/241114si7SKRNAseq/"

# SRR5818141 から SRR5818148 までのファイルをリストする
file_paths <- list.files(path = file_path, pattern = "SRR581814[1-8]_rsem_counts_with_gene_symbols_and_transcript_type.csv", full.names = TRUE)

# ファイルを読み込んでリストに保存
data_list <- lapply(file_paths, read.csv)
# 必要なカラムを取り出すために、すべてのデータをマージ
merged_data <- Reduce(function(x, y) merge(x, y, by="gene_id", all=TRUE), data_list)
# DESeq2パッケージの読み込み
library(DESeq2)
data_combined <- read.csv("C:/Documents/241114si7SKRNAseq/8sample_TPM.csv")
head(data_combined)
# サンプル名を列名として設定
colnames(data_combined)[2:ncol(data_combined)] <- paste0("sample", 1:(ncol(data_combined) - 1))
# 実験群と対照群のサンプル名
control_samples <- c("sample1", "sample2", "sample3", "sample4")
experimental_samples <- c("sample5", "sample6", "sample7", "sample8")

# Fold Change（発現比）の計算
data_combined$fold_change <- rowMeans(data_combined[, experimental_samples]) / rowMeans(data_combined[, control_samples])

# p値の計算 (t検定)
data_combined$p_value <- apply(data_combined[, c(experimental_samples, control_samples)], 1, function(x) t.test(x[1:4], x[5:8])$p.value)

