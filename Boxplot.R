
setwd("C:/Documents/20241228boxplot")
list.dirs("C:/Documents/20241228boxplot")
list.files("C:/Documents/20241228boxplot/csv")
data <- read.csv("C:/Documents/20241228boxplot/csv/si7SK_CDK9_IF.csv")
head(data)
# tidyverse を読み込む
library(tidyverse)

# データを長い形式に変換
data_long <- data %>%
  pivot_longer(cols = c(siCT, si7SK), 
               names_to = "group", 
               values_to = "value")

# group 列を factor として順番を指定 (siCTが左、si7SKが右)
data_long$group <- factor(data_long$group, levels = c("siCT", "si7SK"))

# ggplotでボックスプロットを作成
p <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  labs(x = "Group", y = "Expression", title = "si7SK_CDK9_IF") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 14),  # x軸ラベルの文字サイズ
        axis.text.x = element_text(size = 12),    # x軸目盛り文字サイズ
        axis.title.y = element_text(size = 14),    # y軸ラベルの文字サイズ
        axis.text.y = element_text(size = 12),     # y軸目盛り文字サイズ
        plot.title = element_text(size = 16, face = "bold"))  # タイトルの設定

# グラフを表示
print(p)

# 縦横比7:3で保存
ggsave("C:/Documents/20241228boxplot/boxplot/si7SK_CDK9_IF_boxplot.pdf", plot = p, height = 6, width = 3, units = "in")
###############################################################################################################################
# CSVファイルのディレクトリを指定
csv_dir <- "C:/Documents/20241228boxplot/csv"

# CSVファイルのリストを取得
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

# 各CSVファイルに対して処理を繰り返す
for (csv_file in csv_files) {
  
  # CSVファイルを読み込む
  data <- read.csv(csv_file)
  
  # データを長い形式に変換
  data_long <- data %>%
    pivot_longer(cols = c(siCT, si7SK), 
                 names_to = "group", 
                 values_to = "value")
  
  # group 列を factor として順番を指定 (siCTが左、si7SKが右)
  data_long$group <- factor(data_long$group, levels = c("siCT", "si7SK"))
  
  # ファイル名からプロットタイトルを作成（CSVファイル名を使用）
  plot_title <- gsub(pattern = ".csv$", replacement = "", basename(csv_file))
  
  # ggplotでボックスプロットを作成
  p <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    labs(x = "Group", y = "Expression", title = plot_title) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),  # x軸ラベルの文字サイズ
          axis.text.x = element_text(size = 12),    # x軸目盛り文字サイズ
          axis.title.y = element_text(size = 14),    # y軸ラベルの文字サイズ
          axis.text.y = element_text(size = 12),     # y軸目盛り文字サイズ
          plot.title = element_text(size = 10, face = "bold"))  # タイトルの設定
  
  # グラフを表示
  print(p)
  
  # 保存先ディレクトリとファイル名を設定
  output_file <- paste0("C:/Documents/20241228boxplot/boxplot/", plot_title, "_boxplot.pdf")
  
  # 縦横比6:3で保存
  ggsave(output_file, plot = p, height = 6, width = 3, units = "in")
}
######################################################################################################
#特異点なしﾊﾟﾀﾝ。検定なし
# tidyverse を読み込む
library(tidyverse)

# CSVファイルのディレクトリを指定
csv_dir <- "C:/Documents/20241228boxplot/csv"

# CSVファイルのリストを取得
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

# 各CSVファイルに対して処理を繰り返す
for (csv_file in csv_files) {
  
  # CSVファイルを読み込む
  data <- read.csv(csv_file)
  
  # データを長い形式に変換
  data_long <- data %>%
    pivot_longer(cols = c(siCT, si7SK), 
                 names_to = "group", 
                 values_to = "value")
  
  # group 列を factor として順番を指定 (siCTが左、si7SKが右)
  data_long$group <- factor(data_long$group, levels = c("siCT", "si7SK"))
  
  # ファイル名からプロットタイトルを作成（CSVファイル名を使用）
  plot_title <- gsub(pattern = ".csv$", replacement = "", basename(csv_file))
  
  # ggplotでボックスプロットを作成（特異点を表示しない設定）
  p <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA) +  # 特異点を表示しない
    labs(x = "Group", y = "Expression", title = plot_title) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),  # x軸ラベルの文字サイズ
          axis.text.x = element_text(size = 12),    # x軸目盛り文字サイズ
          axis.title.y = element_text(size = 14),    # y軸ラベルの文字サイズ
          axis.text.y = element_text(size = 12),     # y軸目盛り文字サイズ
          plot.title = element_text(size = 10, face = "bold"))  # タイトルの設定
  
  # グラフを表示
  print(p)
  
  # 保存先ディレクトリとファイル名を設定
  output_file <- paste0("C:/Documents/20241228boxplot/boxplot/", plot_title, "_boxplot_no_outliers.pdf")
  
  # 縦横比6:3で保存
  ggsave(output_file, plot = p, height = 6, width = 3, units = "in")
}


######################################################################################################
# tidyverse を読み込む
library(tidyverse)

# CSVファイルのリストを取得
csv_dir <- "C:/Documents/20241228boxplot/csv"
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

# 各CSVファイルに対して処理を繰り返す
for (csv_file in csv_files) {
  
  # CSVファイルを読み込む
  data <- read.csv(csv_file)
  
  # siCT と si7SK を数値型に変換（文字列として読み込まれている場合に備えて）
  data$siCT <- as.numeric(data$siCT)
  data$si7SK <- as.numeric(data$si7SK)
  
  # データを長い形式に変換
  data_long <- data %>%
    pivot_longer(cols = c(siCT, si7SK), 
                 names_to = "group", 
                 values_to = "value")
  
  # group 列を factor として順番を指定 (siCTが左、si7SKが右)
  data_long$group <- factor(data_long$group, levels = c("siCT", "si7SK"))
  
  # Mann-Whitney U検定（独立した2群の比較）
  # siCT と si7SK の独立した2群を比較する
  wilcox_test_result <- wilcox.test(data$siCT, data$si7SK)
  
  # 検定のp値を取得（非常に小さいp値は"< 0.001"として表示）
  p_value <- ifelse(wilcox_test_result$p.value < 0.001, "< 0.001", round(wilcox_test_result$p.value, 3))
  
  # NAを除外して最大値と最小値を取得
  max_value <- max(data_long$value, na.rm = TRUE)
  min_value <- min(data_long$value, na.rm = TRUE)
  
  # NAがある場合、最大値と最小値が取得できたか確認
  if (is.na(max_value) || is.na(min_value)) {
    print(paste("NA detected in data for file:", csv_file))
    next
  }
  
  print(paste("max:", max_value, "min:", min_value))
  
  # ファイル名からプロットタイトルを作成（CSVファイル名を使用）
  plot_title <- gsub(pattern = ".csv$", replacement = "", basename(csv_file))
  
  # ggplotでボックスプロットを作成
  p <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    labs(x = "Group", y = "Expression", title = plot_title) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),  # x軸ラベルの文字サイズ
          axis.text.x = element_text(size = 12),    # x軸目盛り文字サイズ
          axis.title.y = element_text(size = 14),    # y軸ラベルの文字サイズ
          axis.text.y = element_text(size = 12),     # y軸目盛り文字サイズ
          plot.title = element_text(size = 16, face = "bold")) +  # タイトルの設定
    # 検定結果のp値を表示（位置を動的に決定）
    annotate("text", 
             x = 1.5, 
             y = max_value + (0.1 * (max_value - min_value)), 
             label = paste("p =", p_value), 
             size = 6, 
             hjust = 0.5)
  
  # グラフを表示
  print(p)
  
  # 保存先ディレクトリとファイル名を設定
  output_file <- paste0("C:/Documents/20241228boxplot/boxplot/", plot_title, "_boxplot_with_pvalue.pdf")
  
  # 縦横比6:3で保存
  ggsave(output_file, plot = p, height = 6, width = 3, units = "in")
}
###############################################################################################################
# tidyverse を読み込む
library(tidyverse)

# CSVファイルのリストを取得
csv_dir <- "C:/Documents/20241228boxplot/csv"
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

# 各CSVファイルに対して処理を繰り返す
for (csv_file in csv_files) {
  
  # CSVファイルを読み込む
  data <- read.csv(csv_file)
  
  # siCT と si7SK を数値型に変換（文字列として読み込まれている場合に備えて）
  data$siCT <- as.numeric(data$siCT)
  data$si7SK <- as.numeric(data$si7SK)
  
  # データを長い形式に変換
  data_long <- data %>%
    pivot_longer(cols = c(siCT, si7SK), 
                 names_to = "group", 
                 values_to = "value")
  
  # group 列を factor として順番を指定 (siCTが左、si7SKが右)
  data_long$group <- factor(data_long$group, levels = c("siCT", "si7SK"))
  
  # Mann-Whitney U検定（独立した2群の比較）
  wilcox_test_result <- wilcox.test(data$siCT, data$si7SK)
  
  # 検定のp値を取得（非常に小さいp値は"< 0.001"として表示）
  p_value <- ifelse(wilcox_test_result$p.value < 0.001, "< 0.001", round(wilcox_test_result$p.value, 3))
  
  # siCT と si7SK の n数を計算
  n_siCT <- sum(!is.na(data$siCT))  # siCT の有効な値の数
  n_si7SK <- sum(!is.na(data$si7SK))  # si7SK の有効な値の数
  
  # NAを除外して最大値と最小値を取得
  max_value <- max(data_long$value, na.rm = TRUE)
  min_value <- min(data_long$value, na.rm = TRUE)
  
  # NAがある場合、最大値と最小値が取得できたか確認
  if (is.na(max_value) || is.na(min_value)) {
    print(paste("NA detected in data for file:", csv_file))
    next
  }
  
  # ファイル名からプロットタイトルを作成（CSVファイル名を使用）
  plot_title <- gsub(pattern = ".csv$", replacement = "", basename(csv_file))
  
  # ggplotでボックスプロットを作成
  p <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    labs(x = "Group", y = "Expression", title = plot_title) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),  # x軸ラベルの文字サイズ
          axis.text.x = element_text(size = 12),    # x軸目盛り文字サイズ
          axis.title.y = element_text(size = 14),    # y軸ラベルの文字サイズ
          axis.text.y = element_text(size = 12),     # y軸目盛り文字サイズ
          plot.title = element_text(size = 16, face = "bold")) +  # タイトルの設定
    # 検定結果のp値とn数を右下に表示（位置を動的に決定）
    annotate("text", 
             x = 2.5,  # 右端に近い位置（si7SKの右側）
             y = max_value + 0.1 * (max_value - min_value),  # 最大値より少し上
             label = paste("p =", p_value, 
                           "\nsiCT n =", n_siCT, 
                           "\nsi7SK n =", n_si7SK), 
             size = 3,  # 文字サイズ
             hjust = 1, 
             vjust = 1,
             color = "black")  # 文字の色も指定できます（例: "black")
  
  # グラフを表示
  print(p)
  
  # 保存先ディレクトリとファイル名を設定
  output_file <- paste0("C:/Documents/20241228boxplot/boxplot/", plot_title, "_boxplot_with_pvalue_and_n.pdf")
  
  # 縦横比6:3で保存
  ggsave(output_file, plot = p, height = 6, width = 3, units = "in")
}
####################################################################################################################
# tidyverse を読み込む
library(tidyverse)

# CSVファイルのリストを取得
csv_dir <- "C:/Documents/20241228boxplot/csv_ICE1vsAFF4"
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

# 各CSVファイルに対して処理を繰り返す
for (csv_file in csv_files) {
  
  # CSVファイルを読み込む
  data <- read.csv(csv_file)
  
  # ICE1 と AFF4 を数値型に変換（文字列として読み込まれている場合に備えて）
  data$ICE1 <- as.numeric(data$ICE1)
  data$AFF4 <- as.numeric(data$AFF4)
  
  # データを長い形式に変換
  data_long <- data %>%
    pivot_longer(cols = c(ICE1, AFF4), 
                 names_to = "group", 
                 values_to = "value")
  
  # group 列を factor として順番を指定 (ICE1が左、AFF4が右)
  data_long$group <- factor(data_long$group, levels = c("ICE1", "AFF4"))
  
  # 正規分布を検定（シャピロ・ウィルク検定）
  shapiro_ICE1 <- shapiro.test(data$ICE1)
  shapiro_AFF4 <- shapiro.test(data$AFF4)
  
  # 正規分布でない場合は Mann-Whitney U検定を実行、それ以外は t検定
  if (shapiro_ICE1$p.value < 0.05 | shapiro_AFF4$p.value < 0.05) {
    # Mann-Whitney U検定（独立した2群の比較）
    test_result <- wilcox.test(data$ICE1, data$AFF4)
    test_type <- "Mann-Whitney U test"
    normality_text <- "Non-normal distribution"
  } else {
    # t検定（独立した2群の比較）
    test_result <- t.test(data$ICE1, data$AFF4)
    test_type <- "t-test"
    normality_text <- "Normal distribution"
  }
  
  # 検定のp値を取得（非常に小さいp値は"< 0.001"として表示）
  p_value <- ifelse(test_result$p.value < 0.001, "< 0.001", round(test_result$p.value, 3))
  
  # ICE1 と AFF4 の n数を計算
  n_ICE1 <- sum(!is.na(data$ICE1))  # ICE1 の有効な値の数
  n_AFF4 <- sum(!is.na(data$AFF4))  # AFF4 の有効な値の数
  
  # NAを除外して最大値と最小値を取得
  max_value <- max(data_long$value, na.rm = TRUE)
  min_value <- min(data_long$value, na.rm = TRUE)
  
  # NAがある場合、最大値と最小値が取得できたか確認
  if (is.na(max_value) || is.na(min_value)) {
    print(paste("NA detected in data for file:", csv_file))
    next
  }
  
  # ファイル名からプロットタイトルを作成（CSVファイル名を使用）
  plot_title <- gsub(pattern = ".csv$", replacement = "", basename(csv_file))
  
  # ggplotでボックスプロットを作成
  p <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    labs(x = "Group", y = "Expression", title = plot_title) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),  # x軸ラベルの文字サイズ
          axis.text.x = element_text(size = 12),    # x軸目盛り文字サイズ
          axis.title.y = element_text(size = 14),    # y軸ラベルの文字サイズ
          axis.text.y = element_text(size = 12),     # y軸目盛り文字サイズ
          plot.title = element_text(size = 16, face = "bold")) +  # タイトルの設定
    # 検定結果のp値とn数を右下に表示（位置を動的に決定）
    annotate("text", 
             x = 2.7,  # さらに右に配置
             y = max_value + 0.1 * (max_value - min_value),  # 最大値より少し上
             label = paste(test_type, 
                           "\np =", p_value, 
                           "\nICE1 n =", n_ICE1, 
                           "\nAFF4 n =", n_AFF4,
                           "\nNormality: ", normality_text), 
             size = 3,  # 文字サイズ
             hjust = 1, 
             vjust = 1,
             color = "black")  # 文字の色も指定できます（例: "black")
  
  # グラフを表示
  print(p)
  
  # 保存先ディレクトリとファイル名を設定
  output_file <- paste0("C:/Documents/20241228boxplot/boxplot/", plot_title, "_boxplot_with_pvalue_and_n.pdf")
  
  # 縦横比6:3で保存
  ggsave(output_file, plot = p, height = 6, width = 3, units = "in")
}
#####################################################################################################################
# tidyverse を読み込む
library(tidyverse)

# CSVファイルのリストを取得
csv_dir <- "C:/Documents/20241228boxplot/csv"
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

# 各CSVファイルに対して処理を繰り返す
for (csv_file in csv_files) {
  
  # CSVファイルを読み込む
  data <- read.csv(csv_file)
  
  # siCT と si7SK を数値型に変換（文字列として読み込まれている場合に備えて）
  data$siCT <- as.numeric(data$siCT)
  data$si7SK <- as.numeric(data$si7SK)
  
  # データを長い形式に変換
  data_long <- data %>%
    pivot_longer(cols = c(siCT, si7SK), 
                 names_to = "group", 
                 values_to = "value")
  
  # group 列を factor として順番を指定 (siCTが左、si7SKが右)
  data_long$group <- factor(data_long$group, levels = c("siCT", "si7SK"))
  
  # Mann-Whitney U検定（独立した2群の比較）
  wilcox_test_result <- wilcox.test(data$siCT, data$si7SK)
  
  # 検定のp値を取得（非常に小さいp値は"< 0.001"として表示）
  p_value <- ifelse(wilcox_test_result$p.value < 0.001, "< 0.001", round(wilcox_test_result$p.value, 3))
  
  # siCT と si7SK の n数を計算
  n_siCT <- sum(!is.na(data$siCT))  # siCT の有効な値の数
  n_si7SK <- sum(!is.na(data$si7SK))  # si7SK の有効な値の数
  
  # NAを除外して最大値と最小値を取得
  max_value <- max(data_long$value, na.rm = TRUE)
  min_value <- min(data_long$value, na.rm = TRUE)
  
  # NAがある場合、最大値と最小値が取得できたか確認
  if (is.na(max_value) || is.na(min_value)) {
    print(paste("NA detected in data for file:", csv_file))
    next
  }
  
  # ファイル名からプロットタイトルを作成（CSVファイル名を使用）
  plot_title <- gsub(pattern = ".csv$", replacement = "", basename(csv_file))
  
  # ggplotでボックスプロットを作成
  p <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    labs(x = "Group", y = "Expression", title = plot_title) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),  # x軸ラベルの文字サイズ
          axis.text.x = element_text(size = 12),    # x軸目盛り文字サイズ
          axis.title.y = element_text(size = 14),    # y軸ラベルの文字サイズ
          axis.text.y = element_text(size = 12),     # y軸目盛り文字サイズ
          plot.title = element_text(size = 16, face = "bold")) +  # タイトルの設定
    # 検定結果のp値とn数を右下に表示（位置を動的に決定）
    annotate("text", 
             x = 2.5,  # 右端に近い位置（si7SKの右側）
             y = max_value + 0.1 * (max_value - min_value),  # 最大値より少し上
             label = paste("p =", p_value, 
                           "\nsiCT n =", n_siCT, 
                           "\nsi7SK n =", n_si7SK), 
             size = 3,  # 文字サイズ
             hjust = 1, 
             vjust = 1,
             color = "black")  # 文字の色も指定できます（例: "black")
  
  # グラフを表示
  print(p)
  
  # 保存先ディレクトリとファイル名を設定
  output_file <- paste0("C:/Documents/20241228boxplot/boxplot/", plot_title, "_boxplot_with_pvalue_and_n.pdf")
  
  # 縦横比6:3で保存
  ggsave(output_file, plot = p, height = 6, width = 3, units = "in")
}
####################################################################################################################
# tidyverse を読み込む
library(tidyverse)

# CSVファイルのリストを取得
csv_dir <- "C:/Documents/20250504ELLICE1AFF4"
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

# 各CSVファイルに対して処理を繰り返す
for (csv_file in csv_files) {
  
  # CSVファイルを読み込む
  data <- read.csv(csv_file)
  
  # ICE1 と AFF4 を数値型に変換（文字列として読み込まれている場合に備えて）
  data$ICE1 <- as.numeric(data$ICE1)
  data$AFF4 <- as.numeric(data$AFF4)
  
  # データを長い形式に変換
  data_long <- data %>%
    pivot_longer(cols = c(ICE1, AFF4), 
                 names_to = "group", 
                 values_to = "value")
  
  # group 列を factor として順番を指定 (ICE1が左、AFF4が右)
  data_long$group <- factor(data_long$group, levels = c("ICE1", "AFF4"))
  
  # 正規分布を検定（シャピロ・ウィルク検定）
  shapiro_ICE1 <- shapiro.test(data$ICE1)
  shapiro_AFF4 <- shapiro.test(data$AFF4)
  
  # 正規分布でない場合は Mann-Whitney U検定を実行、それ以外は t検定
  if (shapiro_ICE1$p.value < 0.05 | shapiro_AFF4$p.value < 0.05) {
    # Mann-Whitney U検定（独立した2群の比較）
    test_result <- wilcox.test(data$ICE1, data$AFF4)
    test_type <- "Mann-Whitney U test"
    normality_text <- "Non-normal distribution"
  } else {
    # t検定（独立した2群の比較）
    test_result <- t.test(data$ICE1, data$AFF4)
    test_type <- "t-test"
    normality_text <- "Normal distribution"
  }
  
  # 検定のp値を取得（非常に小さいp値は"< 0.001"として表示）
  p_value <- ifelse(test_result$p.value < 0.001, "< 0.001", round(test_result$p.value, 3))
  
  # ICE1 と AFF4 の n数を計算
  n_ICE1 <- sum(!is.na(data$ICE1))  # ICE1 の有効な値の数
  n_AFF4 <- sum(!is.na(data$AFF4))  # AFF4 の有効な値の数
  
  # NAを除外して最大値と最小値を取得
  max_value <- max(data_long$value, na.rm = TRUE)
  min_value <- min(data_long$value, na.rm = TRUE)
  
  # NAがある場合、最大値と最小値が取得できたか確認
  if (is.na(max_value) || is.na(min_value)) {
    print(paste("NA detected in data for file:", csv_file))
    next
  }
  
  # ファイル名からプロットタイトルを作成（CSVファイル名を使用）
  plot_title <- gsub(pattern = ".csv$", replacement = "", basename(csv_file))
  
  # ggplotでボックスプロットを作成
  p <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    labs(x = "Group", y = "Expression", title = plot_title) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),  # x軸ラベルの文字サイズ
          axis.text.x = element_text(size = 12),    # x軸目盛り文字サイズ
          axis.title.y = element_text(size = 14),    # y軸ラベルの文字サイズ
          axis.text.y = element_text(size = 12),     # y軸目盛り文字サイズ
          plot.title = element_text(size = 16, face = "bold")) +  # タイトルの設定
    # 検定結果のp値とn数を右下に表示（位置を動的に決定）
    annotate("text", 
             x = 2.7,  # さらに右に配置
             y = max_value + 0.1 * (max_value - min_value),  # 最大値より少し上
             label = paste(test_type, 
                           "\np =", p_value, 
                           "\nICE1 n =", n_ICE1, 
                           "\nAFF4 n =", n_AFF4,
                           "\nNormality: ", normality_text), 
             size = 3,  # 文字サイズ
             hjust = 1, 
             vjust = 1,
             color = "black")  # 文字の色も指定できます（例: "black")
  
  # グラフを表示
  print(p)
  
  # 保存先ディレクトリとファイル名を設定
  output_file <- paste0("C:/Documents/20250504ELLICE1AFF4/", plot_title, "_boxplot_with_pvalue_and_n.pdf")
  
  # 縦横比6:3で保存
  ggsave(output_file, plot = p, height = 6, width = 3, units = "in")
}
#############################################################################
C:/Documents/202501117SKkdCOILINcdk9/7skKDcdk9488Coilin555/result_matome_CT4vs7SK3.csv
"E:/molecularbiology/my_Article/supportFigure5/202112027SKsi7SKFISH/matome"
"D:\molecularbiology_20230223\IFData\20231007si7SK_ELL_ENLMED26_NELF\ELL\images\ELL_result\ELL.csv"
# tidyverse を読み込む
library(tidyverse)

# CSVファイルのリストを取得
csv_dir <- "E:/molecularbiology/my_Article/supportFigure5/202112027SKsi7SKFISH/matome"
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

# 各CSVファイルに対して処理を繰り返す
for (csv_file in csv_files) {
  
  # CSVファイルを読み込む
  data <- read.csv(csv_file)
  
  # siCT と si7SK を数値型に変換（文字列として読み込まれている場合に備えて）
  data$siCT <- as.numeric(data$siCT)
  data$si7SK <- as.numeric(data$si7SK)
  
  # データを長い形式に変換
  data_long <- data %>%
    pivot_longer(cols = c(siCT, si7SK), 
                 names_to = "group", 
                 values_to = "value")
  
  # group 列を factor として順番を指定 (siCTが左、si7SKが右)
  data_long$group <- factor(data_long$group, levels = c("siCT", "si7SK"))
  
  # 正規分布を検定（シャピロ・ウィルク検定）
  shapiro_siCT <- shapiro.test(data$siCT)
  shapiro_si7SK <- shapiro.test(data$si7SK)
  
  # 正規分布でない場合は Mann-Whitney U検定を実行、それ以外は t検定
  if (shapiro_siCT$p.value < 0.05 | shapiro_si7SK$p.value < 0.05) {
    # Mann-Whitney U検定（独立した2群の比較）
    test_result <- wilcox.test(data$siCT, data$si7SK)
    test_type <- "Mann-Whitney U test"
  } else {
    # t検定（独立した2群の比較）
    test_result <- t.test(data$siCT, data$si7SK)
    test_type <- "t-test"
  }
  
  # 検定のp値を取得（非常に小さいp値は"< 0.001"として表示）
  p_value <- ifelse(test_result$p.value < 0.001, "< 0.001", round(test_result$p.value, 3))
  
  # siCT と si7SK の n数を計算
  n_siCT <- sum(!is.na(data$siCT))  # siCT の有効な値の数
  n_si7SK <- sum(!is.na(data$si7SK))  # si7SK の有効な値の数
  
  # NAを除外して最大値と最小値を取得
  max_value <- max(data_long$value, na.rm = TRUE)
  min_value <- min(data_long$value, na.rm = TRUE)
  
  # NAがある場合、最大値と最小値が取得できたか確認
  if (is.na(max_value) || is.na(min_value)) {
    print(paste("NA detected in data for file:", csv_file))
    next
  }
  
  # ファイル名からプロットタイトルを作成（CSVファイル名を使用）
  plot_title <- gsub(pattern = ".csv$", replacement = "", basename(csv_file))
  
  # ggplotでボックスプロットを作成
  p <- ggplot(data_long, aes(x = group, y = value, fill = group)) +
    geom_boxplot() +
    labs(x = "Group", y = "Expression", title = plot_title) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),  # x軸ラベルの文字サイズ
          axis.text.x = element_text(size = 12),    # x軸目盛り文字サイズ
          axis.title.y = element_text(size = 14),    # y軸ラベルの文字サイズ
          axis.text.y = element_text(size = 12),     # y軸目盛り文字サイズ
          plot.title = element_text(size = 16, face = "bold")) +  # タイトルの設定
    # 検定結果のp値とn数を右下に表示（位置を動的に決定）
    annotate("text", 
             x = 2.7,  # さらに右に配置
             y = max_value + 0.1 * (max_value - min_value),  # 最大値より少し上
             label = paste(test_type, 
                           "\np =", p_value, 
                           "\nsiCT n =", n_siCT, 
                           "\nsi7SK n =", n_si7SK), 
             size = 3,  # 文字サイズ
             hjust = 1, 
             vjust = 1,
             color = "black")  # 文字の色も指定できます（例: "black")
  
  # グラフを表示
  print(p)
  
  # 保存先ディレクトリとファイル名を設定
  #output_file <- paste0("D:/molecularbiology_20230223/IFData/20231007si7SK_ELL_ENLMED26_NELF/ELL/images/ELL_result/", plot_title, "_boxplot_with_pvalue_and_n.pdf")
  output_file <- paste0(csv_dir, plot_title, "_boxplot_with_pvalue_and_n.pdf")
  
  # 縦横比6:3で保存
  ggsave(output_file, plot = p, height = 6, width = 3, units = "in")
}
############################################################################################################
