
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
    
    # Condition の順番を設定（ICE1 が左、AFF4 が右）
    data$Condition <- factor(data$Condition, levels = c("ICE1", "AFF4"))
    
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
