# 必要なライブラリを読み込み
library(ggplot2)
library(dplyr)
library(readr)

# 対象フォルダのパス
folder_path <- "C:/Documents/qPCR_t_test/Fig3si7SK"

# フォルダ内のすべての CSV ファイルを取得
csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)

# 平均とSEを計算する関数
calculate_summary <- function(data) {
  data %>%
    group_by(Primer, Condition) %>%
    summarise(
      Mean = mean(Value),
      SE = sd(Value) / sqrt(n()),
      .groups = 'drop'
    )
}

# 片側 t 検定を実行する関数
perform_t_test <- function(primer_data) {
  mean_siCT <- mean(primer_data$Value[primer_data$Condition == "siCT"])
  mean_si7SK <- mean(primer_data$Value[primer_data$Condition == "si7SK"])
  
  if (mean_si7SK < mean_siCT) {
    alt_hypothesis <- "greater"
  } else {
    alt_hypothesis <- "less"
  }
  
  p_value <- t.test(Value ~ Condition, data = primer_data, alternative = alt_hypothesis)$p.value
  
  # p値を0.01未満なら "p<0.01" と表示、それ以外は "p=" をつけて表示
  p_value_text <- ifelse(p_value < 0.01, "p<0.01", paste0("p=", format(round(p_value, 4), nsmall = 4)))
  
  return(p_value_text)
}

# 全ての CSV ファイルに対して処理を実行
for (csv_file in csv_files) {
  
  # CSVファイルの読み込み
  data <- read_csv(csv_file)
  
  # Condition の順番を設定（siCT が左、si7SK が右）
  data$Condition <- factor(data$Condition, levels = c("siCT", "si7SK"))
  
  # データの要約を計算
  summary_data <- calculate_summary(data)
  
  # 各 Primer ごとに p 値を計算
  p_values <- data %>%
    group_by(Primer) %>%
    summarise(p_value = perform_t_test(pick(everything())))
  
  # p値をプロット用データに結合
  plot_data <- merge(summary_data, p_values, by = "Primer")
  
  # 全てのPrimerに対して個別にプロットを作成して保存
  unique_primers <- unique(plot_data$Primer)
  
  for (primer_name in unique_primers) {
    
    # 特定のPrimerのデータを抽出
    primer_data <- plot_data[plot_data$Primer == primer_name, ]
    original_data <- data[data$Primer == primer_name, ]
    p_value_text <- primer_data$p_value[1]  # p値を抽出
    
    # ファイル名の元となるCSVファイル名を取得 (拡張子なし)
    base_name <- tools::file_path_sans_ext(basename(csv_file))
    # y軸の最大値を決定（最大値の1.2倍）
    Y_MAX <- max(primer_data$Mean + primer_data$SE) * 1.2
    
    ###### データ点ありバージョン ######
    plot_with_points <- ggplot(primer_data, aes(x = Condition, y = Mean, fill = Condition)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
               color = "black", fill = "black", width = 0.4) +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                    position = position_dodge(width = 0.7), width = 0.2) +
      geom_point(data = original_data, aes(x = Condition, y = Value),
                 position = position_jitter(width = 0.1, height = 0), 
                 size = 5, color = "darkgreen") +
      geom_text(label = p_value_text, 
                x = 1.5, y = Y_MAX * 0.95, size = 10) +  # p値の位置を調整（上限の95%位置）
      scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX)) +  # y軸の上限を固定
      theme_minimal() +
      theme(
        axis.line = element_line(color = "black"),  # 軸を黒くする
        axis.ticks = element_line(color = "black"),# 目盛り線も黒くする
        axis.text.x = element_text(size = 40, color = "black", margin = margin(t = 10)),  # **X軸の文字を黒**
        axis.text.y = element_text(size = 22, color = "black"),  # **Y軸の文字を黒**
        axis.ticks.length = unit(-0.3, "cm")  # **補助目盛りを内向きにする**
      ) +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                    position = position_dodge(width = 0.7), 
                    width = 0.05,  # **T字の横線を短く**
                    size = 0.5) +
      labs(title = paste(base_name, "-", primer_name, "(with points)"),
           x = NULL,  # x軸ラベルを削除
           y = NULL)  # y軸ラベルを削除
    
    
    
    # PDFとして保存 (データ点あり)
    output_file_with_points <- paste0(folder_path, "/", base_name, "_", primer_name, "_with_points.pdf")
    ggsave(filename = output_file_with_points, plot = plot_with_points, width = 6, height = 8)

    ####### データ点なしバージョン ######
    plot_no_points <- ggplot(primer_data, aes(x = Condition, y = Mean, fill = Condition)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
               color = "black", fill = "black", width = 0.4) +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                    position = position_dodge(width = 0.7), width = 0.2) +
      geom_text(label = p_value_text, 
                x = 1.5, y = Y_MAX * 0.95, size = 10) +  # p値の位置を調整（上限の95%位置）
      scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX)) +  # y軸の上限を固定
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 40, color = "black", margin = margin(t = 10)),  # **X軸の文字を黒**
        axis.text.y = element_text(size = 22, color = "black"),  # **Y軸の文字を黒**
        axis.line = element_line(color = "black"),  # 軸を黒くする
        axis.ticks = element_line(color = "black"), # 目盛り線も黒くする
        axis.ticks.length = unit(-0.3, "cm")  # **補助目盛りを内向きにする**
      ) +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                    position = position_dodge(width = 0.7), 
                    width = 0.05,  # **T字の横線を短く**
                    size = 0.5) +  # **線を薄く**
      labs(title = paste(base_name, "-", primer_name, "(no points)"),
           x = NULL,  # x軸ラベルを削除
           y = NULL)  # y軸ラベルを削除
    
    # PDFとして保存 (データ点なし)
    output_file_no_points <- paste0(folder_path, "/", base_name, "_", primer_name, "_no_points.pdf")
    ggsave(filename = output_file_no_points, plot = plot_no_points, width = 6, height = 8)
    
  }
}

# 必要なライブラリを読み込み
library(ggplot2)
library(dplyr)
library(readr)

# 対象の親フォルダのパス
parent_folder_path <- "C:/Documents/qPCR_t_test"

# 親フォルダ内のサブフォルダを取得
sub_folders <- list.dirs(path = parent_folder_path, full.names = TRUE, recursive = FALSE)

# 各サブフォルダに対して処理を実行
for (folder_path in sub_folders) {
  
  # フォルダ内のすべての CSV ファイルを取得
  csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)
  
  # CSVファイルの処理
  for (csv_file in csv_files) {
    
    # CSVファイルの読み込み
    data <- read_csv(csv_file)
    
    # Condition の順番を設定（siCT が左、si7SK が右）
    data$Condition <- factor(data$Condition, levels = c("siCT", "si7SK"))
    
    # データの要約を計算
    summary_data <- calculate_summary(data)
    
    # 各 Primer ごとに p 値を計算
    p_values <- data %>%
      group_by(Primer) %>%
      summarise(p_value = perform_t_test(pick(everything())))
    
    # p値をプロット用データに結合
    plot_data <- merge(summary_data, p_values, by = "Primer")
    
    # 全てのPrimerに対して個別にプロットを作成して保存
    unique_primers <- unique(plot_data$Primer)
    
    for (primer_name in unique_primers) {
      
      # 特定のPrimerのデータを抽出
      primer_data <- plot_data[plot_data$Primer == primer_name, ]
      original_data <- data[data$Primer == primer_name, ]
      p_value_text <- primer_data$p_value[1]  # p値を抽出
      
      # ファイル名の元となるCSVファイル名を取得 (拡張子なし)
      base_name <- tools::file_path_sans_ext(basename(csv_file))
      # y軸の最大値を決定（最大値の1.2倍）
      Y_MAX <- max(primer_data$Mean + primer_data$SE) * 1.2
      
      ###### データ点ありバージョン ######
      plot_with_points <- ggplot(primer_data, aes(x = Condition, y = Mean, fill = Condition)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
                 color = "black", fill = "black", width = 0.3) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                      position = position_dodge(width = 0.7), width = 0.2) +
        geom_point(data = original_data, aes(x = Condition, y = Value),
                   position = position_jitter(width = 0.1, height = 0), 
                   size = 5, color = "darkgreen") +
        geom_text(label = p_value_text, 
                  x = 1.5, y = Y_MAX * 0.95, size = 10, color = "gray40") +  # p値の位置を調整（上限の95%位置）
        scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX)) +  # y軸の上限を固定
        theme_minimal() +
        theme(
          axis.line = element_line(color = "black"),  # 軸を黒くする
          axis.ticks = element_line(color = "black"),# 目盛り線も黒くする
          axis.text.x = element_text(size = 45, color = "black", margin = margin(t = 10)),  # **X軸の文字を黒**
          axis.text.y = element_text(size = 30, color = "black"),  # **Y軸の文字を黒**
          axis.ticks.length = unit(-0.3, "cm")  # **補助目盛りを内向きにする**
        ) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                      position = position_dodge(width = 0.7), 
                      width = 0.05,  # **T字の横線を短く**
                      size = 0.5) +
        labs(title = paste(base_name, "-", primer_name, "(with points)"),
             x = NULL,  # x軸ラベルを削除
             y = NULL)  # y軸ラベルを削除
      
      
      # PDFとして保存 (データ点あり)
      output_file_with_points <- paste0(folder_path, "/", base_name, "_", primer_name, "_with_points.pdf")
      ggsave(filename = output_file_with_points, plot = plot_with_points, width = 6, height = 8)
      
      ####### データ点なしバージョン ######
      plot_no_points <- ggplot(primer_data, aes(x = Condition, y = Mean, fill = Condition)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
                 color = "black", fill = "black", width = 0.3) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                      position = position_dodge(width = 0.7), width = 0.2) +
        geom_text(label = p_value_text, 
                  x = 1.5, y = Y_MAX * 0.95, size = 10, color = "gray40") +  # p値の位置を調整（上限の95%位置）
        scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX)) +  # y軸の上限を固定
        theme_minimal() +
        theme(
          axis.text.x = element_text(size = 45, color = "black", margin = margin(t = 10)),  # **X軸の文字を黒**
          axis.text.y = element_text(size = 30, color = "black"),  # **Y軸の文字を黒**
          axis.line = element_line(color = "black"),  # 軸を黒くする
          axis.ticks = element_line(color = "black"), # 目盛り線も黒くする
          axis.ticks.length = unit(-0.3, "cm")  # **補助目盛りを内向きにする**
        ) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                      position = position_dodge(width = 0.7), 
                      width = 0.05,  # **T字の横線を短く**
                      size = 0.5) +  # **線を薄く**
        labs(title = paste(base_name, "-", primer_name, "(no points)"),
             x = NULL,  # x軸ラベルを削除
             y = NULL)  # y軸ラベルを削除
      
      # PDFとして保存 (データ点なし)
      output_file_no_points <- paste0(folder_path, "/", base_name, "_", primer_name, "_no_points.pdf")
      ggsave(filename = output_file_no_points, plot = plot_no_points, width = 6, height = 8)
      
    }
  }
}
###############################################################################################################
library(tidyverse)

# フォルダのパスを取得
main_folder <- "C:/Documents/qPCR_t_test_pair"  # ここにメインのフォルダパスを指定
sub_folders <- list.dirs(main_folder, recursive = FALSE)  # サブフォルダを取得

# t検定を実行する関数（片側検定）
perform_t_test <- function(sprimer_data) {
  mean_siCT <- mean(primer_data$Value[primer_data$Condition == "siCT"], na.rm = TRUE)
  mean_si7SK <- mean(primer_data$Value[primer_data$Condition == "si7SK"], na.rm = TRUE)
  
  if (mean_si7SK < mean_siCT) {
    alt_hypothesis <- "greater"
  } else {
    alt_hypothesis <- "less"
  }
  
  # 片側 t 検定を実行
  p_value <- t.test(primer_data$Value[primer_data$Condition == "siCT"], 
                    primer_data$Value[primer_data$Condition == "si7SK"], 
                    paired = TRUE, alternative = alt_hypothesis, var.equal = TRUE)$p.value
  
      c
  return(p_value)  # **文字列ではなく数値のまま返す！**
}

# 平均値とSEを計算する関数
calculate_summary <- function(data) {
  data %>%
    group_by(Primer, Condition) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SE = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}

# 各サブフォルダに対して処理を実行
for (folder_path in sub_folders) {
  
  # フォルダ内のすべての CSV ファイルを取得
  csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)
  
  # CSVファイルの処理
  for (csv_file in csv_files) {
    
    # CSVファイルの読み込み
    data <- read_csv(csv_file)
    
    # Condition の順番を設定（siCT が左、si7SK が右）
    data$Condition <- factor(data$Condition, levels = c("siCT", "si7SK"))
    
    # データの要約を計算
    summary_data <- calculate_summary(data)
    
    # 各 Primer ごとに p 値を計算（片側検定）
    p_values <- data %>%
      group_by(Primer) %>%
      summarise(p_value = perform_t_test(pick(everything())))
    
    # p値をプロット用データに結合
    plot_data <- merge(summary_data, p_values, by = "Primer")
    
    # 全てのPrimerに対して個別にプロットを作成して保存
    unique_primers <- unique(plot_data$Primer)
    
    for (primer_name in unique_primers) {
      
      # 特定のPrimerのデータを抽出
      primer_data <- plot_data[plot_data$Primer == primer_name, ]
      original_data <- data[data$Primer == primer_name, ]
      p_value_numeric <- primer_data$p_value[1]  # **数値のまま保持**
      p_value_text <- ifelse(p_value_numeric < 0.01, "***",
                             ifelse(p_value_numeric < 0.05, "**",
                                    ifelse(p_value_numeric < 0.1, "*",
                                           ifelse(p_value_numeric < 0.2, sprintf("p=%.3f", p_value_numeric), "NS"))))# p値を小数点3桁で表示      
      # ファイル名の元となるCSVファイル名を取得 (拡張子なし)
      base_name <- tools::file_path_sans_ext(basename(csv_file))
      
      # y軸の最大値を決定（最大値の1.2倍）
      Y_MAX <- max(primer_data$Mean + primer_data$SE, na.rm = TRUE) * 1.2
      # p値のテキストサイズを条件分岐で設定
      p_value_size <- ifelse(grepl("\\*", p_value_text), 25, 15)  # "*" が含まれていたら25、それ以外なら15
      ###### データ点ありバージョン ######
      plot_with_points <- ggplot(primer_data, aes(x = Condition, y = Mean, fill = Condition)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
                 color = "black", fill = "black", width = 0.3) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                      position = position_dodge(width = 0.7), width = 0.2) +
        geom_point(data = original_data, aes(x = Condition, y = Value),
                   position = position_jitter(width = 0.1, height = 0), 
                   size = 5, color = "darkgreen") +
        geom_text(label = p_value_text, 
                  x = 1.5, y = Y_MAX * 0.95, size = p_value_size, color = "gray40") +  # p値の色をやや濃く
        scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX)) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(size = 45, color = "black", margin = margin(t = 10)),
          axis.text.y = element_text(size = 30, color = "black", margin = margin(r = 10)),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.ticks.length = unit(-0.3, "cm")
        ) +
        labs(title = paste(base_name, "-", primer_name, "(with points)"),
             x = NULL, y = NULL)
      
      # PDFとして保存 (データ点あり)
      output_file_with_points <- paste0(folder_path, "/", base_name, "_", primer_name, "_with_points.pdf")
      ggsave(filename = output_file_with_points, plot = plot_with_points, width = 6, height = 8)
      
      ####### データ点なしバージョン ######
      plot_no_points <- ggplot(primer_data, aes(x = Condition, y = Mean, fill = Condition)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
                 color = "black", fill = "black", width = 0.3) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                      position = position_dodge(width = 0.7), width = 0.2) +
        geom_text(label = p_value_text, 
                  x = 1.5, y = Y_MAX * 0.95, size = p_value_size, color = "gray40") +
        scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX)) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(size = 45, color = "black", margin = margin(t = 10)),
          axis.text.y = element_text(size = 30, color = "black", margin = margin(r = 10)),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.ticks.length = unit(-0.3, "cm")
        ) +
        labs(title = paste(base_name, "-", primer_name, "(no points)"),
             x = NULL, y = NULL)
      
      # PDFとして保存 (データ点なし)
      output_file_no_points <- paste0(folder_path, "/", base_name, "_", primer_name, "_no_points.pdf")
      ggsave(filename = output_file_no_points, plot = plot_no_points, width = 6, height = 8)
      
    }
  }
}
######################################################################################################
# 必要なパッケージを読み込む
library(ggplot2)
library(dplyr)

# データを読み込む
data <- read.csv("C:/Documents/qPCR_t_test/Fig6si7SK/ICE1AFF4_Ser2.csv")

# 平均値とSEを計算
summary_data <- data %>%
  group_by(Primer, Condition) %>%
  summarise(
    Mean = mean(Value, na.rm = TRUE),
    SE = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 片側t検定を実行する関数
perform_t_test <- function(primer_data) {
  mean_siCT <- mean(primer_data$Value[primer_data$Condition == "siCT"], na.rm = TRUE)
  mean_si7SK <- mean(primer_data$Value[primer_data$Condition == "si7SK"], na.rm = TRUE)
  
  alt_hypothesis <- ifelse(mean_si7SK < mean_siCT, "greater", "less")
  
  p_value <- t.test(primer_data$Value[primer_data$Condition == "siCT"], 
                    primer_data$Value[primer_data$Condition == "si7SK"], 
                    paired = FALSE, alternative = alt_hypothesis, var.equal = TRUE)$p.value
  
  return(p_value)
}

# p値を計算
p_values <- data %>%
  group_by(Primer) %>%
  summarise(p_value = perform_t_test(pick(everything())))


# p値のラベルを設定
p_values <- p_values %>%
  mutate(p_label = case_when(
    p_value < 0.01 ~ "***",
    p_value < 0.05 ~ "**",
    p_value < 0.1 ~ "*",
    p_value < 0.2 ~ sprintf("p=%.3f", p_value),
    TRUE ~ "NS"
  ))

# プロット用のデータを作成
plot_data <- merge(summary_data, p_values, by = "Primer")

# `Primer` ごとに個別の PDF を保存
output_dir <- "C:/Documents/qPCR_t_test/Fig6si7SK/plots"
dir.create(output_dir, showWarnings = FALSE)

unique_primers <- unique(plot_data$Primer)

for (primer_name in unique_primers) {
  primer_data <- plot_data[plot_data$Primer == primer_name, ]
  original_data <- data[data$Primer == primer_name, ]
  p_value_text <- primer_data$p_label[1]
  
  # y軸の最大値を決定（余裕を持たせる）
  Y_MAX <- max(primer_data$Mean + primer_data$SE, na.rm = TRUE) * 1.2
  
  # グラフ（データ点あり）
  plot_with_points <- ggplot(primer_data, aes(x = Condition, y = Mean, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
             color = "black", fill = "black", width = 0.3) +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                  position = position_dodge(width = 0.7), width = 0.2) +
    geom_point(data = original_data, aes(x = Condition, y = Value),
               position = position_jitter(width = 0.1, height = 0), 
               size = 3, color = "darkgreen") +
    geom_text(label = p_value_text, 
              x = 1.5, y = Y_MAX * 0.95, size = 8, color = "gray40") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(-0.3, "cm")
    ) +
    labs(title = paste("Primer:", primer_name, "(with points)"),
         x = NULL, y = "Expression Level")
  
  # PDF保存（データ点あり）
  ggsave(filename = file.path(output_dir, paste0(primer_name, "_with_points.pdf")),
         plot = plot_with_points, width = 6, height = 8)
  
  # グラフ（データ点なし）
  plot_no_points <- ggplot(primer_data, aes(x = Condition, y = Mean, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
             color = "black", fill = "black", width = 0.3) +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                  position = position_dodge(width = 0.7), width = 0.2) +
    geom_text(label = p_value_text, 
              x = 1.5, y = Y_MAX * 0.95, size = 10, color = "gray40") +  # p値の位置を調整（上限の95%位置）
    scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX)) +  # y軸の上限を固定
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 45, color = "black", margin = margin(t = 10)),  # **X軸の文字を黒**
      axis.text.y = element_text(size = 30, color = "black"),  # **Y軸の文字を黒**
      axis.line = element_line(color = "black"),  # 軸を黒くする
      axis.ticks = element_line(color = "black"), # 目盛り線も黒くする
      axis.ticks.length = unit(-0.3, "cm")  # **補助目盛りを内向きにする**
    ) +
    geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                  position = position_dodge(width = 0.7), 
                  width = 0.05,  # **T字の横線を短く**
                  size = 0.5) +  # **線を薄く**
    labs(title = paste(base_name, "-", primer_name, "(no points)"),
         x = NULL,  # x軸ラベルを削除
         y = NULL)  # y軸ラベルを削除
  
  # PDFとして保存 (データ点なし)
  output_file_no_points <- paste0(folder_path, "/", base_name, "_", primer_name, "_no_points.pdf")
  ggsave(filename = output_file_no_points, plot = plot_no_points, width = 6, height = 8)
  # PDF保存（データ点なし）
  ggsave(filename = file.path(output_dir, paste0(primer_name, "_no_points.pdf")),
         plot = plot_no_points, width = 6, height = 8)
}

print("PDFファイルが保存されました。")


plot_no_points <- ggplot(primer_data, aes(x = Condition, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
           color = "black", fill = "black", width = 0.3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(label = p_value_text, 
            x = 1.5, y = Y_MAX * 0.95, size = 10, color = "gray40") +  # p値の位置を調整（上限の95%位置）
  scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX)) +  # y軸の上限を固定
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 45, color = "black", margin = margin(t = 10)),  # **X軸の文字を黒**
    axis.text.y = element_text(size = 30, color = "black"),  # **Y軸の文字を黒**
    axis.line = element_line(color = "black"),  # 軸を黒くする
    axis.ticks = element_line(color = "black"), # 目盛り線も黒くする
    axis.ticks.length = unit(-0.3, "cm")  # **補助目盛りを内向きにする**
  ) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                position = position_dodge(width = 0.7), 
                width = 0.05,  # **T字の横線を短く**
                size = 0.5) +  # **線を薄く**
  labs(title = paste(base_name, "-", primer_name, "(no points)"),
       x = NULL,  # x軸ラベルを削除
       y = NULL)  # y軸ラベルを削除

# PDFとして保存 (データ点なし)
output_file_no_points <- paste0(folder_path, "/", base_name, "_", primer_name, "_no_points.pdf")
ggsave(filename = output_file_no_points, plot = plot_no_points, width = 6, height = 8)
########################################################################################

library(ggplot2)
library(dplyr)
library(readr)
# 対象の親フォルダのパス
parent_folder_path <- "C:/Documents/tsuika_bar"

# 親フォルダ内のサブフォルダを取得
sub_folders <- list.dirs(path = parent_folder_path, full.names = TRUE, recursive = FALSE)

calculate_summary <- function(df) {
  df %>%
    group_by(Primer, Condition) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SE = sd(Value, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
}

perform_t_test <- function(df) {
  tryCatch({
    test_result <- t.test(Value ~ Condition, data = df, var.equal = TRUE)  # Studentのt検定
    return(test_result$p.value)
  }, error = function(e) {
    return(NA)
  })
}


# 各サブフォルダに対して処理を実行
for (folder_path in sub_folders) {
  
  # フォルダ内のすべての CSV ファイルを取得
  csv_files <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)
  
  # CSVファイルの処理
  for (csv_file in csv_files) {
    
    # CSVファイルの読み込み
    data <- read_csv(csv_file)
    
    # Condition の順番を設定（siCT が左、si7SK が右）
    data$Condition <- factor(data$Condition, levels = c("siCT", "si7SK"))
    
    # データの要約を計算
    summary_data <- calculate_summary(data)
    
    # 各 Primer ごとに p 値を計算
    p_values <- data %>%
      group_by(Primer) %>%
      summarise(p_value = perform_t_test(pick(everything())))
    
    # p値に基づく記号または値の割り当て
    p_values <- p_values %>%
      mutate(significance = ifelse(p_value < 0.01, "***",
                                   ifelse(p_value < 0.05, "**",
                                          ifelse(p_value < 0.1, "*",
                                                 ifelse(p_value < 0.2, sprintf("p=%.3f", p_value), "NS")))))
    
    # p値をプロット用データに結合
    plot_data <- merge(summary_data, p_values, by = "Primer")
    
    # 全てのPrimerに対して個別にプロットを作成して保存
    unique_primers <- unique(plot_data$Primer)
    
    for (primer_name in unique_primers) {
      
      # 特定のPrimerのデータを抽出
      primer_data <- plot_data[plot_data$Primer == primer_name, ]
      original_data <- data[data$Primer == primer_name, ]
      p_value_text <- primer_data$significance[1]  # significanceを抽出
      
      # ファイル名の元となるCSVファイル名を取得 (拡張子なし)
      base_name <- tools::file_path_sans_ext(basename(csv_file))
      # y軸の最大値を決定（最大値の1.2倍）
      Y_MAX <- max(primer_data$Mean + primer_data$SE) * 1.2
      
      ###### データ点ありバージョン ######
      plot_with_points <- ggplot(primer_data, aes(x = Condition, y = Mean, fill = Condition)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
                 color = "black", fill = "black", width = 0.3) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                      position = position_dodge(width = 0.7), width = 0.2) +
        geom_point(data = original_data, aes(x = Condition, y = Value),
                   position = position_jitter(width = 0.1, height = 0), 
                   size = 5, color = "darkgreen") +
        geom_text(label = p_value_text, 
                  x = 1.5, y = Y_MAX * 0.95, size = 10, color = "gray40") +  # p値の位置を調整（上限の95%位置）
        scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX)) +  # y軸の上限を固定
        theme_minimal() +
        theme(
          axis.line = element_line(color = "black"),  # 軸を黒くする
          axis.ticks = element_line(color = "black"),# 目盛り線も黒くする
          axis.text.x = element_text(size = 45, color = "black", margin = margin(t = 10)),  # **X軸の文字を黒**
          axis.text.y = element_text(size = 30, color = "black"),  # **Y軸の文字を黒**
          axis.ticks.length = unit(-0.3, "cm")  # **補助目盛りを内向きにする**
        ) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                      position = position_dodge(width = 0.7), 
                      width = 0.05,  # **T字の横線を短く**
                      size = 0.5) +
        labs(title = paste(base_name, "-", primer_name, "(with points)"),
             x = NULL,  # x軸ラベルを削除
             y = NULL)  # y軸ラベルを削除
      
      
      # PDFとして保存 (データ点あり)
      output_file_with_points <- paste0(folder_path, "/", base_name, "_", primer_name, "_with_points.pdf")
      ggsave(filename = output_file_with_points, plot = plot_with_points, width = 6, height = 8)
      p_value_size <- ifelse(grepl("\\*", p_value_text), 25, 15)
      ####### データ点なしバージョン ######
      plot_no_points <- ggplot(primer_data, aes(x = Condition, y = Mean, fill = Condition)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
                 color = "black", fill = "black", width = 0.3) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                      position = position_dodge(width = 0.7), width = 0.2) +
        geom_text(label = p_value_text, 
                  x = 1.5, y = Y_MAX * 0.95, size = p_value_size, color = "gray40") +
        scale_y_continuous(expand = c(0, 0), limits = c(0, Y_MAX)) +  # y軸の上限を固定
        theme_minimal() +
        theme(
          axis.text.x = element_text(size = 45, color = "black", margin = margin(t = 10)),  # **X軸の文字を黒**
          axis.text.y = element_text(size = 30, color = "black"),  # **Y軸の文字を黒**
          axis.line = element_line(color = "black"),  # 軸を黒くする
          axis.ticks = element_line(color = "black"), # 目盛り線も黒くする
          axis.ticks.length = unit(-0.3, "cm")  # **補助目盛りを内向きにする**
        ) +
        geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE),
                      position = position_dodge(width = 0.7), 
                      width = 0.05,  # **T字の横線を短く**
                      size = 0.5) +  # **線を薄く**
        labs(title = paste(base_name, "-", primer_name, "(no points)"),
             x = NULL,  # x軸ラベルを削除
             y = NULL)  # y軸ラベルを削除
      
      # PDFとして保存 (データ点なし)
      output_file_no_points <- paste0(folder_path, "/", base_name, "_", primer_name, "_no_points.pdf")
      ggsave(filename = output_file_no_points, plot = plot_no_points, width = 6, height = 8)
      
    }
  }
}
print(paste("Saving PDF to:", output_file_with_points))
ggsave(filename = output_file_with_points, plot = plot_with_points, width = 6, height = 8)
# この直前に
base_name <- tools::file_path_sans_ext(basename(csv_file))
print(paste("Base name:", base_name))  # これを追加

output_file_with_points <- paste0(folder_path, "/", base_name, "_", primer_name, "_with_points.pdf")
print(paste("Saving PDF to:", output_file_with_points))  # これを追加


