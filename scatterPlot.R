library(biomaRt)
setwd("C:/Documents/20241119Ser2Pol2scatterplot")
# 1. BioMart に接続
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # 人の場合

# 2. 遺伝子名リストの読み込み
gene_list <- read.csv("C:/Documents/20241119Ser2Pol2scatterplot/gene_list.csv")  # ファイル名を指定
gene_names <- gene_list$gene_name  # "gene_name" 列から遺伝子名を取得

# 3. 取得したい属性を指定
attributes <- c(
  "ensembl_gene_id",      # 遺伝子の安定ID (ENSG...)
  "external_gene_name",   # 遺伝子名
  "gene_biotype",         # 遺伝子の種類 (例: protein_coding, pseudogene)
  "chromosome_name",      # 染色体
  "start_position",       # 遺伝子の開始位置
  "end_position"          # 遺伝子の終了位置
)

# 4. 遺伝子情報を取得
results <- getBM(
  attributes = attributes,
  filters = "external_gene_name",  # フィルターとして遺伝子名を指定
  values = gene_names,             # 遺伝子名リストを渡す
  mart = ensembl
)

# 5. 結果を確認
print(head(results))

# 6. CSV に保存
write.csv(results, "gene_information_with_types.csv", row.names = FALSE)
  packageVersion("pheatmap")

library(biomaRt)# 1. BioMart に接続
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # 人の場合 (Homo sapiens)

# 2. 取得したい遺伝子名のリスト
gene_list <- c("LOC124904835", "MIR520D", "RNVU1-19", "RNU1-1")

# 3. 取得する属性（columns の指定）
attributes <- c(
  "ensembl_gene_id",  # Gene stable ID
  "external_gene_name",  # 遺伝子名
  "gene_biotype",  # 遺伝子の種類 (Gene Biotype)
  "chromosome_name",  # 染色体名
  "start_position",  # 遺伝子開始位置
  "end_position"  # 遺伝子終了位置
)

# 4. 遺伝子情報を取得
results <- getBM(
  attributes = attributes,
  filters = "external_gene_name",  # フィルターを遺伝子名で指定
  values = gene_list,  # フィルターに渡す値
  mart = ensembl
)

# 5. 結果の確認
print(results)

# 必要であれば CSV に保存
write.csv(results, "gene_information.csv", row.names = FALSE)

# 遺伝子名に基づいて情報を取得（ensemble_gene_idも含む）
gene_info <- getBM(attributes = c("external_gene_name", "gene_biotype", "ensembl_gene_id"),
                   filters = "external_gene_name", 
                   values = gene_names, 
                   mart = ensembl)

~/NGS/gencode.v45.annotation.gtf


# 結果表示
print(gene_info)

library(biomaRt)
# EnsemblのBioMartサーバーに接続
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # ヒトの場合

# 必要な遺伝子情報を取得
genes <- getBM(
  attributes = c("chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id", "gene_biotype"),
  mart = ensembl
)

# データの確認
head(genes)
# "chr"を付ける
genes$chromosome_name <- paste0("chr", genes$chromosome_name)

# 開始位置を0-basedに変換
genes$start_position <- genes$start_position - 1

# BED形式に整形（必要なカラムのみ選択）
bed_data <- genes[, c("chromosome_name", "start_position", "end_position", "ensembl_gene_id")]

# 小数点以下を削除（必要に応じて）
bed_data$ensembl_gene_id <- sub("\\..*", "", bed_data$ensembl_gene_id)

# BEDファイルとして保存
write.table(
  bed_data,
  file = "C:/Documents/20241119Ser2Pol2scatterplot/all_genes.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)
# ENSG順に並べ替え
bed_data <- bed_data[order(bed_data$ensembl_gene_id), ]

# BEDファイルとして保存
write.table(
  bed_data,
  file = "C:/Documents/20241119Ser2Pol2scatterplot/all_genes.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# biomaRt パッケージを読み込む
library(biomaRt)

# CSVファイルを読み込む
df <- read.csv("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/SRR11440180.trim_output.csv")

# CSVファイルから遺伝子名（ENSG）列を抽出
gene_ids <- df$ENSG

# Ensembl Mart への接続（https://www.ensembl.org）
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "https://www.ensembl.org")

# データセットの選択（遺伝子情報用）
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# 取得したい情報（遺伝子名、遺伝子タイプ）
attributes <- c("ensembl_gene_id", "external_gene_name", "gene_biotype")

# ENSG をキーにして gene_name と gene_biotype を一括取得
gene_info <- getBM(attributes = attributes, 
                   filters = "ensembl_gene_id", 
                   values = gene_ids, 
                   mart = ensembl)

# 取得した遺伝子情報を確認
head(gene_info)

# 元のデータフレームに gene_type をマージ
df_with_type <- merge(df, gene_info, by.x = "ENSG", by.y = "ensembl_gene_id", all.x = TRUE)

# 結果を新しいCSVとして保存
write.csv(df_with_type, "C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/SRR11440180.trim_output_with_gene_type.csv", row.names = FALSE)

# biomaRt パッケージを読み込む
library(biomaRt)

# ファイルリストを指定
files <- c("SRR11440180.trim_output.csv", 
           "SRR11440181.trim_output.csv", 
           "SRR11440186.trim_output.csv", 
           "SRR11440187.trim_output.csv", 
           "SRR11440190.trim_output.csv", 
           "SRR11440191.trim_output.csv", 
           "SRR11440194.trim_output.csv", 
           "SRR11440195.trim_output.csv")

# Ensembl Mart への接続（https://www.ensembl.org）
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "https://www.ensembl.org")

# データセットの選択（遺伝子情報用）
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# 取得したい情報（遺伝子名、遺伝子タイプ）
attributes <- c("ensembl_gene_id", "external_gene_name", "gene_biotype")

# すべてのファイルに対して処理を繰り返す
for (file in files) {
  # ファイルのフルパスを指定
  file_path <- paste0("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/", file)
  
  # CSVファイルを読み込む
  df <- read.csv(file_path)
  
  # CSVファイルから遺伝子名（ENSG）列を抽出
  gene_ids <- df$ENSG
  
  # ENSG をキーにして gene_name と gene_biotype を一括取得
  gene_info <- getBM(attributes = attributes, 
                     filters = "ensembl_gene_id", 
                     values = gene_ids, 
                     mart = ensembl)
  
  # 元のデータフレームに gene_name と gene_biotype をマージ
  df_with_type <- merge(df, gene_info, by.x = "ENSG", by.y = "ensembl_gene_id", all.x = TRUE)
  
  # 新しいCSVファイルのパスを作成
  output_file <- paste0("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/", gsub(".csv", "_with_gene_type.csv", file))
  
  # 結果を新しいCSVとして保存
  write.csv(df_with_type, output_file, row.names = FALSE)
  
  # 終了メッセージ（進捗確認用）
  cat("Processed file:", file, "\n")
}

# 必要なライブラリ
library(dplyr)

# ファイルリストを指定
files <- c("SRR11440180.trim_output_with_gene_type.csv", 
           "SRR11440181.trim_output_with_gene_type.csv", 
           "SRR11440186.trim_output_with_gene_type.csv", 
           "SRR11440187.trim_output_with_gene_type.csv", 
           "SRR11440190.trim_output_with_gene_type.csv", 
           "SRR11440191.trim_output_with_gene_type.csv", 
           "SRR11440194.trim_output_with_gene_type.csv", 
           "SRR11440195.trim_output_with_gene_type.csv")

# ベースディレクトリを指定
base_dir <- "C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/"

# 最終的なデータフレームを格納するリスト
mean_columns <- list()

# 各ファイルから Mean 列を抽出
for (file in files) {
  # ファイルのフルパスを作成
  file_path <- paste0(base_dir, file)
  
  # ファイルを読み込む
  df <- read.csv(file_path)
  
  # ファイル名から SRR 番号を取得
  srr_number <- sub(".*(SRR[0-9]+).*", "\\1", file)
  
  # Mean 列をリストに追加（列名として SRR 番号を設定）
  mean_columns[[srr_number]] <- df$Mean
}

# リストをデータフレームに変換して列を並べる
result_df <- as.data.frame(mean_columns)

# 結果を確認

head(result_df)
# 必要に応じて保存
write.csv(result_df, "C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/Mean_combined_sorted.csv", row.names = FALSE)
# SRR11440180.trim_output_with_gene_type.csvを読み込む
annotated_file <- paste0(base_dir, "SRR11440180.trim_output_with_gene_type.csv")
annotated_data <- read.csv(annotated_file)

# 必要な3列 (ENSG, external_gene_name, gene_biotype) を抽出
annotation <- annotated_data %>%
  select(ENSG, external_gene_name, gene_biotype)

# result_df に annotation を結合（行数が一致していることを前提）
final_df <- cbind(annotation, result_df)

# 結果を確認
print(head(final_df))
# 列名の変更
colnames(final_df)[4:11] <- c(
  "SRR11440180_input1", 
  "SRR11440181_input2", 
  "SRR11440186_Pol2_1", 
  "SRR11440187_Pol2_2", 
  "SRR11440190_Ser2_1", 
  "SRR11440191_Ser2_2", 
  "SRR11440194_Ser5_1", 
  "SRR11440195_Ser5_2"
)

# CSVとして保存
write.csv(final_df, "C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/Mean_combined_renamed.csv", row.names = FALSE)
final_df <- read.csv("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/Mean_combined_renamed.csv")
head(final_df)


# 平均値列を計算
final_df <- final_df %>%
  mutate(
    input_Ave = rowMeans(select(., SRR11440180_input1, SRR11440181_input2), na.rm = TRUE),
    Pol2_Ave = rowMeans(select(., SRR11440186_Pol2_1, SRR11440187_Pol2_2), na.rm = TRUE),
    Ser2_Ave = rowMeans(select(., SRR11440190_Ser2_1, SRR11440191_Ser2_2), na.rm = TRUE),
    Ser5_Ave = rowMeans(select(., SRR11440194_Ser5_1, SRR11440195_Ser5_2), na.rm = TRUE)
  )

# 必要な列のみを選択
summary_df <- final_df %>%
  select(ENSG, external_gene_name, gene_biotype, input_Ave, Pol2_Ave, Ser2_Ave, Ser5_Ave)

# 新しいデータフレームをCSVとして保存
write.csv(summary_df, "C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/Mean_summary_with_averages.csv", row.names = FALSE)

########################################################################################################
##ここまでで、データフレームが出来上がった。UMAPの書き方がしたから始まる
final_df <- read.csv("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/Mean_summary_with_averages.csv")#ここまでのまとめ
library(umap)
library(ggplot2)

# 最初に、遺伝子発現データ部分を抽出
expression_data <- final_df[, 4:ncol(final_df)]  # 発現データ（8サンプル分）

# UMAPを実行
umap_result <- umap(expression_data)

# 結果をデータフレームに変換
umap_df <- data.frame(UMAP1 = umap_result$layout[, 1], 
                      UMAP2 = umap_result$layout[, 2])

# サンプル名を追加
umap_df$Sample <- colnames(expression_data)

# UMAPプロットを作成
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Sample)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "UMAP of 8 Samples", x = "UMAP1", y = "UMAP2") +
  theme(legend.title = element_blank())

# 必要なライブラリのインストールと読み込み
if (!require(umap)) install.packages("umap", dependencies = TRUE)
if (!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)

library(umap)
library(ggplot2)

# 最初に、遺伝子発現データ部分を抽出
expression_data <- final_df[, 4:ncol(final_df)]  # 発現データ（8サンプル分）
head(expression_data)
# サンプル名を手動で設定
sample_names <- c("input_1", "input_2", "Pol2_1", "Pol2_2", "Ser2_1", "Ser2_2", "Ser5_1", "Ser5_2")

# UMAPを実行
umap_result <- umap(expression_data)

# 結果をデータフレームに変換
umap_df <- data.frame(UMAP1 = umap_result$layout[, 1], 
                      UMAP2 = umap_result$layout[, 2])
head(umap_result)
# サンプル名を追加
umap_df$Sample <- sample_names

# サンプル名をumap_dfに追加する
# サンプル名の長さがumap_resultの行数と一致することを確認する
if (length(sample_names) == nrow(umap_df)) {
  umap_df$Sample <- sample_names
} else {
  stop("サンプル名の数とUMAP結果の行数が一致しません。")
}


# UMAPプロットを作成
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Sample)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "UMAP of 8 Samples", x = "UMAP1", y = "UMAP2") +
  theme(legend.title = element_blank())



####################################################################################################
df <- read.csv("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/Mean_summary_with_averages.csv")#ここまでのまとめ
setwd("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed")
library(dplyr)
library(ggplot2)
head(df)
# 基本の散布図
ggplot(df, aes(x = Pol2_Ave, y = Ser2_Ave)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Scatter Plot: Pol2_Ave vs Ser2_Ave",
    x = "Pol2_Ave",
    y = "Ser2_Ave"
  )
# 平均列の分散を計算
summary_df <- df %>%
  mutate(variance = rowSums(select(., input_Ave, Pol2_Ave, Ser2_Ave, Ser5_Ave)^2) / 4)

# 分散の上位10%を選択
high_variance_genes <- summary_df %>%
  filter(variance > quantile(variance, 0.9))
head(high_variance_genes)
#３次元マッピング
install.packages("plotly")
library(plotly)

######## 3次元散布図を作成###################
library(plotly)
library(rgl)
install.packages("rgl")
head(high_variance_genes)

# 3D 散布図を作成
plot <- plot_ly(
  data = high_variance_genes,
  x = ~Pol2_Ave,
  y = ~Ser2_Ave,
  z = ~Ser5_Ave,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 5, color = ~variance, colorscale = 'Viridis', showscale = TRUE)
)

install.packages("webshot")
webshot::install_phantomjs()
install.packages("htmlwidgets")
library(plotly)

# 3D 散布図を作成
plot <- plot_ly(
  data = high_variance_genes,
  x = ~Pol2_Ave,
  y = ~Ser2_Ave,
  z = ~Ser5_Ave,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 5, color = ~variance, colorscale = 'Viridis', showscale = TRUE)
)
# SVG ファイルとして保存
htmlwidgets::saveWidget(plot, "scatter3d_plot.html")
system("python -m plotly convert scatter3d_plot.html --to svg") 
# 必要なパッケージをロード
library(png)
library(grid)

# PNG 画像を読み込み
img <- readPNG("3DSer2.png")

# PDF ファイルに出力
pdf("scatter3d_plot.pdf")

# PDF に PNG 画像を描画
grid.raster(img)

# PDF 出力を完了
dev.off()


#########################################################################################
library(plotly)

# 回帰モデルの作成
model <- lm(Ser5_Ave ~ Pol2_Ave + Ser2_Ave, data = high_variance_genes)

# グリッドデータを作成
x_vals <- seq(min(high_variance_genes$Pol2_Ave), max(high_variance_genes$Pol2_Ave), length.out = 30)
y_vals <- seq(min(high_variance_genes$Ser2_Ave), max(high_variance_genes$Ser2_Ave), length.out = 30)
grid <- expand.grid(Pol2_Ave = x_vals, Ser2_Ave = y_vals)

# 回帰平面の z 値を予測
grid$Ser5_Ave <- predict(model, newdata = grid)

# 3D散布図
plot <- plot_ly(
  data = high_variance_genes,
  x = ~Pol2_Ave,
  y = ~Ser2_Ave,
  z = ~Ser5_Ave,
  type = 'scatter3d',
  mode = 'markers',
  marker = list(size = 5, color = ~variance, colorscale = 'Viridis', showscale = TRUE)
)

# 回帰平面を追加

plot <- plot %>%
  add_surface(
    x = matrix(grid$Pol2_Ave, nrow = length(x_vals), ncol = length(y_vals)),
    y = matrix(grid$Ser2_Ave, nrow = length(x_vals), ncol = length(y_vals)),
    z = matrix(grid$Ser5_Ave, nrow = length(x_vals), ncol = length(y_vals)),
    opacity = 0.5,  # 透明度を設定
    colorscale = list(c(0, 1), c("lightgray", "lightgray"))  # 薄いグレー
  )


# プロットを表示
plot

##################################################################################################

# 選択した遺伝子の散布図
ggplot(high_variance_genes, aes(x = Pol2_Ave, y = Ser2_Ave)) +
  geom_point(alpha = 0.6, color = "grey") +
  theme_minimal() +
  labs(
    title = "High-Variance Genes: Pol2_Ave vs Ser2_Ave",
    x = "Pol2_Ave",
    y = "Ser2_Ave"
  )

# 分散上位10%に対する散布図を作成
ggplot(high_variance_df, aes(x = Pol2_Ave, y = Ser2_Ave)) +
  # 全体の点を灰色で描画
  geom_point(alpha = 0.6, color = "gray") +
  # snRNA の点を赤色で上書き
  geom_point(data = subset(high_variance_df, gene_biotype == "snRNA"), 
             aes(x = Pol2_Ave, y = Ser2_Ave), 
             color = "red", 
             alpha = 0.8, 
             size = 1.5) +
  # 軸ラベルとテーマの設定
  theme_minimal() +
  labs(
    title = "Scatter Plot: High Variance Genes (Highlight snRNA)",
    x = "Pol2_Ave",
    y = "Ser2_Ave"
  )

# 選択した遺伝子を `high_variance_genes` に追加
selected_genes <- selected_genes_info$ENSG
highlighted_genes <- high_variance_genes %>%
  filter(ENSG %in% selected_genes) %>%
  mutate(highlight = "selected")  # 選択した遺伝子には "selected" タグを付ける

# 他の遺伝子には "other" タグを付ける
high_variance_genes <- high_variance_genes %>%
  filter(!(ENSG %in% selected_genes)) %>%
  mutate(highlight = "other")

# 高分散遺伝子全体を結合
combined_data <- bind_rows(highlighted_genes, high_variance_genes)

# 散布図の作成
ggplot(combined_data, aes(x = Pol2_Ave, y = Ser2_Ave)) +
  geom_point(aes(color = highlight, size = highlight), alpha = 0.6) +
  scale_color_manual(values = c("selected" = "red", "other" = "grey")) + # 色を設定
  scale_size_manual(values = c("selected" = 2, "other" = 1)) +  # サイズを設定
  theme_minimal() +
  labs(
    title = "High-Variance Genes: Pol2_Ave vs Ser2_Ave",
    x = "Pol2_Ave",
    y = "Ser2_Ave"
  ) +
  theme(legend.position = "none") # 凡例を非表示に



############################################################################################
library(ggplot2)
head(df)
# 基本の散布図
ggplot(df, aes(x = Pol2_Ave, y = Ser2_Ave)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Scatter Plot: Pol2_Ave vs Ser2_Ave",
    x = "Pol2_Ave",
    y = "Ser2_Ave"
  )
# 平均列の分散を計算
summary_df <- df %>%
  mutate(variance = rowSums(select(., input_Ave, Pol2_Ave, Ser2_Ave, Ser5_Ave)^2) / 4)

# 分散の上位10%を選択
high_variance_genes <- summary_df %>%
  filter(variance > quantile(variance, 0.9))

# 選択した遺伝子の散布図
ggplot(high_variance_genes, aes(x = Pol2_Ave, y = Ser2_Ave)) +
  geom_point(alpha = 0.6, color = "grey") +
  theme_minimal() +
  labs(
    title = "High-Variance Genes: Pol2_Ave vs Ser2_Ave",
    x = "Pol2_Ave",
    y = "Ser2_Ave"
  )

# 分散上位10%に対する散布図を作成
ggplot(high_variance_df, aes(x = Pol2_Ave, y = Ser2_Ave)) +
  # 全体の点を灰色で描画
  geom_point(alpha = 0.6, color = "gray") +
  # snRNA の点を赤色で上書き
  geom_point(data = subset(high_variance_df, gene_biotype == "snRNA"), 
             aes(x = Pol2_Ave, y = Ser2_Ave), 
             color = "red", 
             alpha = 0.8, 
             size = 1.5) +
  # 軸ラベルとテーマの設定
  theme_minimal() +
  labs(
    title = "Scatter Plot: High Variance Genes (Highlight snRNA)",
    x = "Pol2_Ave",
    y = "Ser2_Ave"
  )

# 選択した遺伝子を `high_variance_genes` に追加
selected_genes <- selected_genes_info$ENSG
highlighted_genes <- high_variance_genes %>%
  filter(ENSG %in% selected_genes) %>%
  mutate(highlight = "selected")  # 選択した遺伝子には "selected" タグを付ける

# 他の遺伝子には "other" タグを付ける
high_variance_genes <- high_variance_genes %>%
  filter(!(ENSG %in% selected_genes)) %>%
  mutate(highlight = "other")

# 高分散遺伝子全体を結合
combined_data <- bind_rows(highlighted_genes, high_variance_genes)

# 散布図の作成
ggplot(combined_data, aes(x = Pol2_Ave, y = Ser2_Ave)) +
  geom_point(aes(color = highlight, size = highlight), alpha = 0.6) +
  scale_color_manual(values = c("selected" = "red", "other" = "grey")) + # 色を設定
  scale_size_manual(values = c("selected" = 2, "other" = 1)) +  # サイズを設定
  theme_minimal() +
  labs(
    title = "High-Variance Genes: Pol2_Ave vs Ser2_Ave",
    x = "Pol2_Ave",
    y = "Ser2_Ave"
  ) +
  theme(legend.position = "none") # 凡例を非表示に



############################################################################################
library(ggplot2)
head(df)
# 基本の散布図
ggplot(df, aes(x = Pol2_Ave, y = Ser5_Ave)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Scatter Plot: Pol2_Ave vs Ser5_Ave",
    x = "Pol2_Ave",
    y = "Ser5_Ave"
  )
# 平均列の分散を計算
summary_df <- df %>%
  mutate(variance = rowSums(select(., input_Ave, Pol2_Ave, Ser2_Ave, Ser5_Ave)^2) / 4)

# 分散の上位10%を選択
high_variance_genes <- summary_df %>%
  filter(variance > quantile(variance, 0.9))

# 選択した遺伝子の散布図
ggplot(high_variance_genes, aes(x = Pol2_Ave, y = Ser5_Ave)) +
  geom_point(alpha = 0.6, color = "grey") +
  theme_minimal() +
  labs(
    title = "High-Variance Genes: Pol2_Ave vs Ser5_Ave",
    x = "Pol2_Ave",
    y = "Ser5_Ave"
  )

# 分散上位10%に対する散布図を作成
ggplot(high_variance_df, aes(x = Pol2_Ave, y = Ser2_Ave)) +
  # 全体の点を灰色で描画
  geom_point(alpha = 0.6, color = "gray") +
  # snRNA の点を赤色で上書き
  geom_point(data = subset(high_variance_df, gene_biotype == "snRNA"), 
             aes(x = Pol2_Ave, y = Ser2_Ave), 
             color = "red", 
             alpha = 0.8, 
             size = 1.5) +
  # 軸ラベルとテーマの設定
  theme_minimal() +
  labs(
    title = "Scatter Plot: High Variance Genes (Highlight snRNA)",
    x = "Pol2_Ave",
    y = "Ser2_Ave"
  )

# 選択した遺伝子を `high_variance_genes` に追加
selected_genes <- selected_genes_info$ENSG
highlighted_genes <- high_variance_genes %>%
  filter(ENSG %in% selected_genes) %>%
  mutate(highlight = "selected")  # 選択した遺伝子には "selected" タグを付ける

# 他の遺伝子には "other" タグを付ける
high_variance_genes <- high_variance_genes %>%
  filter(!(ENSG %in% selected_genes)) %>%
  mutate(highlight = "other")

# 高分散遺伝子全体を結合
combined_data <- bind_rows(highlighted_genes, high_variance_genes)

# 散布図の作成
ggplot(combined_data, aes(x = Pol2_Ave, y = Ser5_Ave)) +
  geom_point(aes(color = highlight, size = highlight), alpha = 0.6) +
  scale_color_manual(values = c("selected" = "red", "other" = "grey")) + # 色を設定
  scale_size_manual(values = c("selected" = 2, "other" = 1)) +  # サイズを設定
  theme_minimal() +
  labs(
    title = "High-Variance Genes: Pol2_Ave vs Ser5_Ave",
    x = "Pol2_Ave",
    y = "Ser5_Ave"
  ) +
  theme(legend.position = "none") # 凡例を非表示に
################################################################################################
# 必要なライブラリを読み込み
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
valid_types <- c("lincRNA", "miRNA", "misc_RNA", "ribozyme", "rRNA", "scaRNA", "snoRNA", "snRNA", "sRNA", "TEC", "vault_RNA")

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

# UMAPの結果がV1, V2でない場合、列名を変更する
colnames(umap_data) <- c("UMAP1", "UMAP2")

# UMAPの結果に遺伝子タイプの情報を追加
umap_data$group <- top_genes$group

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
                                "TEC" = "grey", 
                                "vault_RNA" = "darkgreen")) +
  theme_minimal() +
  labs(title = "UMAP of top 10% variance genes", x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "right")

head(df)

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
head(selected_genes_info)
############################################################################################

# データフレームの読み込み（仮定）
df_1 <- read.csv("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/Mean_summary_with_averages.csv")
setwd("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed")
# 分散計算: 各サンプルの変動（std）を計算
variance <- apply(df_1[, c("input_Ave", "Pol2_Ave", "Ser2_Ave")], 1, var)

# 分散上位10％の遺伝子を抽出
threshold <- quantile(variance, 0.9)
top_genes <- df_1[variance >= threshold, ]
top_genes <- high_variance_genes

# UMAPの適用
umap_model <- umap(top_genes[, c("input_Ave", "Pol2_Ave", "Ser5_Ave" )])

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
# UMAPの結果がV1, V2でない場合、列名を変更する
colnames(umap_data) <- c("UMAP1", "UMAP2")

# UMAPの結果に遺伝子タイプの情報を追加
umap_data$group <- top_genes$group

# UMAPプロット
ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point(alpha = 1, size = 1.5) +  # 点を小さめにして視認性を確保
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

head(df_1)


umap_data$ENSG <- top_genes$ENSG# UMAP結果に遺伝子ID（ENSG）を追加
selected_genes <- selected_genes_info$ENSG# 選択した遺伝子のUMAP座標を抽出
selected_umap_data <- umap_data[umap_data$ENSG %in% selected_genes, ]# selected_genesに対応するtop_genesのインデックスを取得
head(selected_umap_data)# 選択した遺伝子のUMAP座標を確認

# UMAPのプロットを作成し、選択した遺伝子を強調表示
ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = group)) +
  geom_point(alpha = 0.5, size = 1.5) +  # 通常のプロット
  geom_point(data = selected_umap_data, aes(x = UMAP1, y = UMAP2), 
             color = "red", size = 3, alpha = 1) +  # 選択した遺伝子を赤で大きく表示
  scale_color_manual(values = c("others" = "grey", 
                                "lincRNA" = "blue", 
                                "miRNA" = "green", 
                                "misc_RNA" = "purple", 
                                "ribozyme" = "orange", 
                                "rRNA" = "pink", 
                                "scaRNA" = "brown", 
                                "snoRNA" = "yellow", 
                                "snRNA" = "darkred", 
                                "sRNA" = "magenta", 
                                "TEC" = "black",
                                "protein_coding"="cyan",
                                "vault_RNA" = "darkgreen")) +
  theme_minimal() +
  labs(title = "UMAP of top 10% variance genes", x = "UMAP1", y = "UMAP2") +
  theme(legend.position = "right")
###################################################################################
head(selected_genes_info)
# ggplot2パッケージを読み込む
library(ggplot2)

# gene_biotypeのカテゴリごとにカウントする
biotype_counts <- table(selected_genes_info$gene_biotype)

# カウントデータをデータフレームに変換
biotype_df <- as.data.frame(biotype_counts)
colnames(biotype_df) <- c("gene_biotype", "count")

# 円グラフを描画
ggplot(biotype_df, aes(x = "", y = count, fill = gene_biotype)) +
  geom_bar(stat = "identity", width = 1) + 
  coord_polar(theta = "y") +  # 円グラフにする
  theme_void() +  # グラフの背景を無くす
  labs(title = "Gene Biotype Distribution") +  # タイトルを追加
  theme(legend.title = element_blank())  # 凡例のタイトルを非表示

# gene_biotypeのカテゴリごとにカウントする
biotype_counts <- table(selected_genes_info$gene_biotype)

# カウントデータをデータフレームに変換
biotype_df <- as.data.frame(biotype_counts)
colnames(biotype_df) <- c("gene_biotype", "count")

# 合計数を計算
biotype_df$fraction <- biotype_df$count / sum(biotype_df$count)

# 割合が多い順に並べ替え
biotype_df <- biotype_df[order(biotype_df$fraction, decreasing = TRUE), ]

# gene_biotypeの順番を逆に指定
biotype_df$gene_biotype <- factor(biotype_df$gene_biotype, levels = c("miRNA", "lncRNA", "snoRNA", "protein_coding", "snRNA"))

# 色の順番と変更
ggplot(biotype_df, aes(x = "", y = fraction, fill = gene_biotype, label = gene_biotype)) +
  geom_bar(stat = "identity", width = 1, color = "black") +  # 枠線を追加
  coord_polar(theta = "y") +  # 円グラフにする
  theme_void() +  # 背景を無くす
  labs(title = "Gene Biotype Distribution") +  # タイトル
  theme(legend.title = element_blank()) +  # 凡例のタイトルを非表示
  scale_fill_manual(values = c("miRNA" = "darkgreen", 
                               "lncRNA" = "lightgreen", 
                               "snoRNA" = "lightblue", 
                               "protein_coding" = "mediumblue", 
                               "snRNA" = "darkblue")) +  # 色指定（順番通りに）
  geom_text(aes(label = paste0(gene_biotype, ": ", scales::percent(fraction))),
            position = position_stack(vjust = 0.5),  # ラベルを円の中に配置
            color = ifelse(biotype_df$gene_biotype == "snRNA" | biotype_df$gene_biotype == "protein_coding"| biotype_df$gene_biotype == "miRNA", "white", "black"))  # 白文字はsnRNAとsnoRNA

R.version

################################################################################################
COILChAP <- read.csv("20221021_CoilinChAPpeakCall_Annotated.csv")
head(COILChAP)
head(selected_genes_info)
install.packages("ggvenn")  # ggvennパッケージのインストール
install.packages("venn")    # vennパッケージのインストール
# COILChAPのENSG情報
ensg_coil_chap <- unique(COILChAP$feature)
# selected_genes_infoのENSG情報
ensg_pol2umapselect <- unique(selected_genes_info$ENSG)

# それぞれのサンプル数
coil_chap_size <- length(ensg_coil_chap)
pol2umapselect_size <- length(ensg_pol2umapselect)

# ggvennを使ってベン図を作成
library(ggvenn)

# データをリスト形式に変換
venn_data <- list(
  COILChAP = ensg_coil_chap,
  Pol2UMAPselect = ensg_pol2umapselect
)

library(ggvenn)

# データをリスト形式に変換
venn_data <- list(
  COILChAP = ensg_coil_chap,
  Pol2UMAPselect = ensg_pol2umapselect
head(venn_data)

# ggvennでベン図を作成
ggvenn(venn_data, 
       set_name_size = 10, 
       stroke_size = 0.5,
       fill_color = c("skyblue", "orange"))

install.packages("VennDiagram")



# 遺伝子数を設定
coil_chap_genes <- 110  # COILChAPの遺伝子数
pol2umapselect_genes <- 31  # Pol2UMAPselectの遺伝子数
intersection_genes <- 48  # 交差部分の遺伝子数

# ベン図用のデータ
venn_data <- list(
  COILChAP = coil_chap_genes,
  Pol2UMAPselect = pol2umapselect_genes
)

# Venn図の作成
venn.plot <- venn.diagram(
  x = venn_data,  # データセット
  category.names = c("COILChAP", "Pol2UMAPselect"),  # ラベル名
  filename = NULL,  # ファイルに保存しない場合はNULL
  fill = c("skyblue", "orange"),  # 円の色
  alpha = 0.5,  # 円の透明度
  cex = 2,  # 数字ラベルのフォントサイズ
  cat.cex = 2,  # カテゴリラベルのフォントサイズ
  cat.pos = c(0, 180),  # カテゴリラベルの位置
  cat.dist = 0.05,  # カテゴリラベルと円との距離
  lwd = 2,  # 線の太さ
  fontface = "bold",  # ラベルのフォントスタイル
  fontfamily = "sans",  # フォントのファミリs
  main = "Venn Diagram of Gene Overlap"  # タイトル
)

# 結果をプロット
grid.draw(venn.plot)
##########################################
library(VennDiagram)
# COILChAP の ENSG 列を抽出
coil_genes <- COILChAP$feature  # "feature" 列に ENSG 番号が含まれている場合
# selected_genes_info の ENSG 列を抽出
selected_genes <- selected_genes_info$ENSG  # "ENSG" 列

# NAが含まれていないか確認
summary(coil_genes_set)
summary(selected_genes_set)

# NAが含まれていれば除去
coil_genes_set <- na.omit(coil_genes_set)
selected_genes_set <- na.omit(selected_genes_set)
# ベン図用のデータセット
venn_data <- list(
  COILChAP = coil_genes_set,  # COILChAPの遺伝子セット

  Pol2UMAPselect = selected_genes_set  # selected_genes_infoの遺伝子セット
)

# Venn図の作成
venn.diagram(
  x = venn_data,  # 遺伝子セット
  filename = "ida_r_venn_diagram_01.svg",  # 保存するファイル名
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
##################################################################################
#peakアノテーション
install.packages("GenomicRanges")
install.packages("ChIPseeker")
BiocManager::install("ChIPseeker")
a
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("ChIPseeker")
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# BEDファイルを読み込む
data.bed <- read.table("C:/Documents/20241119Ser2Pol2scatterplot/SRX595703.10.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(data.bed) <- c("chr", "start", "end", "peak_name", "score", "strand", 
                        "signal_value", "p_value", "q_value", "peak_id")
# GenomicRangesオブジェクトを作成
peak_gr <- GRanges(seqnames = data.bed$chr,
                   ranges = IRanges(start = data.bed$start, end = data.bed$end))
# TxDbオブジェクトを使って遺伝子アノテーション
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# ピークを遺伝子領域にアノテート
peak_annotation <- annotatePeak(peak_gr, TxDb = txdb, tssRegion = c(-300, 300))
class(peak_annotation)
# csAnnoオブジェクトをデータフレームに変換
peak_annotation_df <- as.data.frame(peak_annotation)
head(peak_annotation_df)

# biomaRtがインストールされていない場合、インストールします
# install.packages("biomaRt")

# biomaRtをロード
library(biomaRt)

# Ensemblマートを設定
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# transcriptId から ENSG ID と gene_name を取得
gene_info <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "ensembl_transcript_id"),
                   filters="ensembl_transcript_id",
                   values=peak_annotation_df$transcriptId,
                   mart=ensembl)

# 結果を確認
head(gene_info)
# peak_annotation_df の先頭行を確認
head(peak_annotation_df$transcriptId)
# transcriptId の先頭行を確認
head(peak_annotation_df$transcriptId)

# biomaRt で取得する際に必要な形式が合っているか確認
# biomaRt クエリを再実行して、確認
gene_info <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "ensembl_transcript_id"),
                   filters="ensembl_transcript_id",  # filtersにensembl_transcript_idを指定
                   values=peak_annotation_df$transcriptId,
                   mart=ensembl)

# 結果を確認
head(gene_info)

  library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# transcriptId から遺伝子ID（ENSG）を取得
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_info_txdb <- select(txdb, keys=peak_annotation_df$transcriptId, columns=c("GENEID", "TXNAME"), keytype="TXID")

# 結果を確認
head(gene_info_txdb)
# TxDb で利用できるキーを確認
keys(txdb, keytype="TXID")
library(biomaRt)

# Ensemblのマートを使用
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# transcriptId から ENSG ID と遺伝子名を取得
gene_info <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "ensembl_transcript_id"),
                   filters="ensembl_transcript_id",
                   values=peak_annotation_df$transcriptId,
                   mart=ensembl)

# 結果を確認
head(gene_info)

