setwd("C:/Documents/20250504SER2ICE1AFF4")
install.packages("fs")  # 最初だけ必要
library(fs)

dir_tree("C:/Documents/20250509IFlineScan_Fig3_Fig4/")
y <- read.csv("C:/Documents/20250509IFlineScan_Fig3_Fig4/Ser2_matome/Values_7SK_Coilin.csv")
head(y)


library("tidyverse")
library("ggplot2")
library("ggsci")


g1 <- ggplot(y, aes(x = Distance_.microns.)) +
  geom_line(aes(y = COILIN, color = "COILIN"), size = 1.5) +
  geom_line(aes(y = SER2, color = "SER2"), size = 1.5) +
  scale_color_manual(values = c("COILIN" = "#F8766D", "SER2" = "#00BFC4")) +
  labs(x = "Distance (microns)", y = "Signal Intensity", color = "Signal Type") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "#F5F5F5", color = NA),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "#E0E0E0", linetype = "dotted"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "#F5F5F5", color = NA)
  )  # ここで縦横比1:2（縦:横）を設定

# プロット表示
plot(g1)


# PDFに保存（縦4 x 横7のサイズ）
ggsave("Ser2_matome_Values_7SK_Coilin_line_scan_plot.pdf", plot = g1, width = 7, height = 3, units = "in")
########################################################
# CSVファイルパス
csv_path <- "C:/Documents/20250509IFlineScan_Fig3_Fig4/Ser2_matome/Values_7SK_Coilin.csv"

# データ読み込み
y <- read.csv(csv_path)

# ライブラリ
library(tidyverse)
library(ggsci)

# 3列目の列名（何でも対応できるように取得）
signal_col <- names(y)[3]

# ロング形式に変換：2列目（COILIN）と3列目（可変）を対象に
y_long <- y %>%
  pivot_longer(cols = c("COILIN", all_of(signal_col)),
               names_to = "SignalType", values_to = "Intensity")

# 色を自動でつける（固定色 + 自動色でもよい）
line_colors <- c("COILIN" = "#F8766D", setNames("#00BFC4", signal_col))

# プロット
g1 <- ggplot(y_long, aes(x = Distance_.microns., y = Intensity, color = SignalType)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = line_colors) +
  labs(x = "Distance (microns)", y = "Intensity", color = "Signal Type") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "#F5F5F5", color = NA),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "#E0E0E0", linetype = "dotted"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "#F5F5F5", color = NA)
  )


# 保存ファイル名作成（自動）
folder_path <- dirname(csv_path)
folder_name <- basename(folder_path)
csv_file_name <- tools::file_path_sans_ext(basename(csv_path))
pdf_file_name <- paste0(folder_name, "_", csv_file_name, "_line_scan_plot.pdf")
pdf_path <- file.path(folder_path, pdf_file_name)

# PDF保存
ggsave(pdf_path, plot = g1, width = 7, height = 3, units = "in")
############################################################
# 必要なライブラリ
library(tidyverse)
library(ggsci)

# 親ディレクトリ（ここを基準にして探索）
parent_dir <- "C:/Documents/20250509IFlineScan_Fig3_Fig4/"

# 2階層下のディレクトリ内の CSV ファイルを探索
csv_files <- list.files(parent_dir, recursive = TRUE, full.names = TRUE, pattern = "Values_.*\\.csv")

# CSV ファイルを1つずつ処理
for (csv_path in csv_files) {
  
  # データ読み込み
  y <- read.csv(csv_path)
  
  # 3列目の列名（何でも対応できるように取得）
  signal_col <- names(y)[3]
  
  # COILIN 列が存在するか確認し、存在する場合に処理を行う
  if ("COILIN" %in% names(y)) {
    # ロング形式に変換：COILIN と 3列目（可変）を対象に
    y_long <- y %>%
      pivot_longer(cols = c("COILIN", all_of(signal_col)),
                   names_to = "SignalType", values_to = "Intensity")
    
    # 色を自動でつける（固定色 + 自動色でもよい）
    line_colors <- c("COILIN" = "#F8766D", setNames("#00BFC4", signal_col))
    
    # プロット
    g1 <- ggplot(y_long, aes(x = Distance_.microns., y = Intensity, color = SignalType)) +
      geom_line(size = 1.5) +
      scale_color_manual(values = line_colors) +
      labs(x = "Distance (microns)", y = "Intensity", color = "Signal Type") +
      theme_minimal(base_size = 14) +
      theme(
        panel.background = element_rect(fill = "#F5F5F5", color = NA),
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "#E0E0E0", linetype = "dotted"),
        legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(fill = "#F5F5F5", color = NA)
      )
    
    # 保存ファイル名作成（自動）
    folder_path <- dirname(csv_path)
    folder_name <- basename(folder_path)
    csv_file_name <- tools::file_path_sans_ext(basename(csv_path))
    pdf_file_name <- paste0(folder_name, "_", csv_file_name, "_line_scan_plot.pdf")
    pdf_path <- file.path(folder_path, pdf_file_name)
    
    # PDF保存
    ggsave(pdf_path, plot = g1, width = 7, height = 3, units = "in")
  } else {
    warning(paste("COILIN 列が見つかりません: ", csv_path))
  }
}

