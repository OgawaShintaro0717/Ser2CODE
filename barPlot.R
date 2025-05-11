"C:/Documents/20250118ICEAFFpol2ChIPseq/ICE1_AFF4_RTrato.csv"

data <- read.csv("C:/Documents/20241228boxplot/csv_ICE1vsAFF4_bar/ICE1vsAFF4_RNAseq_Log2(RTratio).csv")
head(data)
"C:/Documents/20241228boxplot/boxplot/"
# 必要なライブラリを読み込み
library(tidyverse)

# CSVデータを読み込み
data <- read.csv("C:/Documents/20250118ICEAFFpol2ChIPseq/ICE1_AFF4_RTrato.csv")
# データを長い形式に変換
data_long <- data %>%
  pivot_longer(cols = c(ICE1, AFF4), 
               names_to = "group", 
               values_to = "Log2_RTratio")

# 各行を一意に識別するIDを作成して、ペアごとにグループを指定
data_long <- data_long %>%
  mutate(pair_id = rep(1:nrow(data), each = 2))  # 1行ごとにペアを作成

# ggplotで直線を作成：ICE1とAFF4をペアで結ぶ
p <- ggplot(data_long, aes(x = group, y = Log2_RTratio, group = pair_id)) +
  geom_line(aes(color = pair_id), linewidth = 1) +  # 直線を描画
  geom_point(aes(color = pair_id), size = 2) +  # 各点をプロット
  labs(x = "Group", y = "Log2(RTratio)", title = "ICE1 vs AFF4") +
  scale_x_discrete(limits = c("ICE1", "AFF4")) +  # x軸をICE1, AFF4の順に変更
  theme_minimal() +
  theme(axis.title.x = element_text(size = 14),  # x軸ラベルの文字サイズ
        axis.text.x = element_text(size = 12),    # x軸目盛り文字サイズ
        axis.title.y = element_text(size = 14),    # y軸ラベルの文字サイズ
        axis.text.y = element_text(size = 12),     # y軸目盛り文字サイズ
        plot.title = element_text(size = 16, face = "bold"))  # タイトルの設定

# グラフを表示
print(p)

# 保存先ディレクトリとファイル名を設定
output_file_2 <- "C:/Documents/20250118ICEAFFpol2ChIPseq/ICE1_AFF4_RTrato_LinePlot_2.pdf"

# グラフを保存
ggsave(output_file_2, plot = p, height = 6, width = 4, units = "in")
######################################################################################

# 必要なライブラリを読み込み
library(tidyverse)

# CSVファイルのディレクトリを指定
csv_dir <- "C:/Documents/20250319_RTratio"

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
                 values_to = "Log2_RTratio")
  
  # 各行を一意に識別するIDを作成して、ペアごとにグループを指定
  data_long <- data_long %>%
    mutate(pair_id = rep(1:nrow(data), each = 2))  # 1行ごとにペアを作成
  
  # ggplotで直線を作成：siCTとsi7SKをペアで結ぶ
  p <- ggplot(data_long, aes(x = group, y = Log2_RTratio, group = pair_id)) +
    geom_line(aes(color = pair_id), linewidth = 1) +  # 直線を描画
    geom_point(aes(color = pair_id), size = 2) +  # 各点をプロット
    labs(x = "Group", y = "Log2(RTratio)", title = basename(csv_file)) +
    scale_x_discrete(limits = c("siCT", "si7SK")) +  # x軸をsiCT, si7SKの順に変更
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),  # x軸ラベルの文字サイズ
          axis.text.x = element_text(size = 12),    # x軸目盛り文字サイズ
          axis.title.y = element_text(size = 14),    # y軸ラベルの文字サイズ
          axis.text.y = element_text(size = 12),     # y軸目盛り文字サイズ
          plot.title = element_text(size = 16, face = "bold"))  # タイトルの設定
  
  # グラフを表示
  print(p)
  
  # 保存先ディレクトリとファイル名を設定
  plot_title <- gsub(pattern = ".csv$", replacement = "_plot.pdf", basename(csv_file))
  output_file <- file.path("C:/Documents/20241228boxplot/boxplot", plot_title)
  
  # 縦横比6:4で保存
  ggsave(output_file, plot = p, height = 6, width = 4, units = "in")
}
##########################################################################

# 必要なライブラリを読み込み
library(tidyverse)

# CSVファイルのディレクトリを指定
csv_dir <- "C:/Documents/20250319_RTratio"

# CSVファイルのリストを取得
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

# 各CSVファイルに対して処理を繰り返す
for (csv_file in csv_files) {
  
  # CSVファイルを読み込む
  data <- read.csv(csv_file)
  
  # データを長い形式に変換
  data_long <- data %>%
    pivot_longer(cols = c(NCM460, HCT116), 
                 names_to = "group", 
                 values_to = "Log2_RTratio")
  
  # 各行を一意に識別するIDを作成して、ペアごとにグループを指定
  data_long <- data_long %>%
    mutate(pair_id = rep(1:nrow(data), each = 2))  # 1行ごとにペアを作成
  
  # ggplotで直線を作成：NCM460とHCT116をペアで結ぶ
  p <- ggplot(data_long, aes(x = group, y = Log2_RTratio, group = pair_id)) +
    geom_line(aes(color = pair_id), linewidth = 1) +  # 直線を描画
    geom_point(aes(color = pair_id), size = 2) +  # 各点をプロット
    labs(x = "Group", y = "Log2(RTratio)", title = basename(csv_file)) +
    scale_x_discrete(limits = c("NCM460", "HCT116")) +  # x軸をNCM460, HCT116の順に変更
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),  # x軸ラベルの文字サイズ
          axis.text.x = element_text(size = 12),    # x軸目盛り文字サイズ
          axis.title.y = element_text(size = 14),    # y軸ラベルの文字サイズ
          axis.text.y = element_text(size = 12),     # y軸目盛り文字サイズ
          plot.title = element_text(size = 16, face = "bold"))  # タイトルの設定
  
  # グラフを表示
  print(p)
  
  # 保存先ディレクトリとファイル名を設定
  plot_title <- gsub(pattern = ".csv$", replacement = "_plot.pdf", basename(csv_file))
  output_file <- file.path("C:/Documents/20241228boxplot/boxplot", plot_title)
  
  # 縦横比6:4で保存
  ggsave(output_file, plot = p, height = 6, width = 4, units = "in")
}
