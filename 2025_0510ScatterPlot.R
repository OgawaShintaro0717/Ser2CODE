df <- read.csv("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/Mean_summary_with_averages.csv")#ここまでのまとめ
df_1 <- read.csv("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/SRR11440180.trim_output_with_gene_type.csv")
df_2 <- read.csv("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/Mean_summary_with_averages_size.csv")#ここまでのまとめ

setwd("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed")
library(dplyr)
library(ggplot2)
head(df)
head(df_1)
dim(df)
dim(df_1)
dim(df_2)
 使用する列を指定
value_cols <- c("input_Ave", "Pol2_Ave", "Ser2_Ave", "Ser5_Ave")

# 分散計算と抽出
df_highvar <- df %>%
  mutate(variance = rowMeans(across(all_of(value_cols), ~ .^2))) %>%
  filter(variance > quantile(variance, 0.9))

# 分散計算と抽出
df_highvar <- df_2 %>%
  mutate(variance = rowMeans(across(all_of(value_cols), ~ .^2))) %>%
  filter(variance > quantile(variance, 0.9))

head(df_highvar)
dim(df_highvar)
df_3 <- read.csv("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/selected_genes_info.csv")
head(df_3)
dim(df_3)
df_highvar <- df_highvar %>%
  mutate(is_in_df3 = external_gene_name %in% df_3$external_gene_name)


library(ggplot2)



write.csv(df_highvar, "C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/df_highvar.csv", row.names = FALSE)

df_highvar_no_df3 <- df_highvar_no_df3 %>% select(-is_in_df3)
head(df_highvar_no_df3)

df_3 <- df_3 %>% select(-external_gene_name)
df_3 <- df_3 %>% select(-gene_biotype)
head(df_3)
df_highvar_no_df3$is_in_df3 <- df_highvar_no_df3$ENSG %in% df_3$ENSG
sum(df_highvar_no_df3$is_in_df3)

######################################################################################
#上記DFのまとめ
write.csv(df_highvar_no_df3, "C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/df_highvar_size.csv", row.names = FALSE)
df_highvar_no_df3 <- read.csv("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/df_highvar_size.csv")

ggplot(df_highvar_no_df3, aes(x = Pol2_Ave, y = Ser2_Ave)) +
  geom_point(aes(color = is_in_df3, size = 1/Size), alpha = 1) +  # Sizeを逆数にして点を調整
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#2C3E50")) +  # TRUEは赤、FALSEは青
  scale_size_continuous(range = c(1, 10)) +  # サイズを適切な範囲に設定
  labs(
    x = "Pol2 Average",
    y = "log2(Ser2 Average)",
    title = "Pol2 vs log2(Ser2) Scatter Plot (High Variance Genes)",
    color = "In df3",
    size = "Gene Size"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0.5)
  )

ggplot(df_highvar_no_df3, aes(x = Pol2_Ave, y = log2(Ser2_Ave))) +
  geom_point(aes(color = is_in_df3), size = 4, alpha = 0.8) +  # 点の大きさを4に設定
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#2C3E50")) +  # TRUEは赤、FALSEは青
  labs(
    x = "Pol2 Average",
    y = "log2(Ser2 Average)",  # 縦軸ラベルを変更
    title = "Pol2 vs log2(Ser2) Scatter Plot (High Variance Genes)",
    color = "In df3"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0.5)
  )
#################################################################################################
# 回帰直線の計算
lm_red <- lm(log2(Ser2_Ave) ~ Pol2_Ave, data = subset(df_highvar_no_df3, is_in_df3 == TRUE))
lm_black <- lm(log2(Ser2_Ave) ~ Pol2_Ave, data = subset(df_highvar_no_df3, is_in_df3 == FALSE))

# 直線の式を文字列で作成
eq_red <- paste("y = ", round(coef(lm_red)[2], 2), "x + ", round(coef(lm_red)[1], 2), sep = "")
eq_black <- paste("y = ", round(coef(lm_black)[2], 2), "x + ", round(coef(lm_black)[1], 2), sep = "")

# グラフ作成
ggplot(df_highvar_no_df3, aes(x = Pol2_Ave, y = log2(Ser2_Ave))) +
  geom_point(aes(color = is_in_df3), size = 2, alpha = 1) +  # 点の大きさを4に設定
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#2C3E50")) +  # TRUEは赤、FALSEは青
  geom_smooth(data = subset(df_highvar_no_df3, is_in_df3 == TRUE), 
              aes(x = Pol2_Ave, y = log2(Ser2_Ave)), 
              method = "lm", color = "red", se = FALSE) +  # 赤の回帰直線
  geom_smooth(data = subset(df_highvar_no_df3, is_in_df3 == FALSE), 
              aes(x = Pol2_Ave, y = log2(Ser2_Ave)), 
              method = "lm", color = "black", se = FALSE) +  # 黒の回帰直線
  geom_text(aes(x = max(Pol2_Ave), y = log2(Ser2_Ave)[which.max(Pol2_Ave)], 
                label = eq_red), 
            color = "red", hjust = 1, vjust = 1) +  # 赤直線の式
  geom_text(aes(x = max(Pol2_Ave), y = log2(Ser2_Ave)[which.max(Pol2_Ave)], 
                label = eq_black), 
            color = "black", hjust = 1, vjust = 0) +  # 黒直線の式
  labs(
    x = "Pol2 Average",
    y = "log2(Ser2 Average)",  # 縦軸ラベルを変更
    title = "Pol2 vs log2(Ser2) Scatter Plot (High Variance Genes)",
    color = "In df3"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0.5)
  )
# PDFとして保存
ggsave("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/scatterplot.pdf", width = 8, height = 6)

# 回帰直線の計算
lm_red <- lm(Ser2_Ave ~ Pol2_Ave, data = subset(df_highvar_no_df3, is_in_df3 == TRUE))
lm_black <- lm(Ser2_Ave ~ Pol2_Ave, data = subset(df_highvar_no_df3, is_in_df3 == FALSE))

# 直線の式を文字列で作成
eq_red <- paste("y = ", round(coef(lm_red)[2], 2), "x + ", round(coef(lm_red)[1], 2), sep = "")
eq_black <- paste("y = ", round(coef(lm_black)[2], 2), "x + ", round(coef(lm_black)[1], 2), sep = "")

# グラフ作成（y軸をSer2_Aveに変更）
ggplot(df_highvar_no_df3, aes(x = Pol2_Ave, y = Ser2_Ave)) +
  geom_point(aes(color = is_in_df3), size = 2, alpha = 1) +  # 点の大きさを4に設定
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#2C3E50")) +  # TRUEは赤、FALSEは青
  geom_smooth(data = subset(df_highvar_no_df3, is_in_df3 == TRUE), 
              aes(x = Pol2_Ave, y = Ser2_Ave), 
              method = "lm", color = "red", se = FALSE) +  # 赤の回帰直線
  geom_smooth(data = subset(df_highvar_no_df3, is_in_df3 == FALSE), 
              aes(x = Pol2_Ave, y = Ser2_Ave), 
              method = "lm", color = "black", se = FALSE) +  # 黒の回帰直線
  geom_text(aes(x = max(Pol2_Ave), y = Ser2_Ave[which.max(Pol2_Ave)], 
                label = eq_red), 
            color = "red", hjust = 1, vjust = 1) +  # 赤直線の式
  geom_text(aes(x = max(Pol2_Ave), y = Ser2_Ave[which.max(Pol2_Ave)], 
                label = eq_black), 
            color = "black", hjust = 1, vjust = 0) +  # 黒直線の式
  labs(
    x = "Pol2 Average",
    y = "Ser2 Average",  # 縦軸ラベルを変更
    title = "Pol2 vs Ser2 Scatter Plot (High Variance Genes)",
    color = "In df3"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0.5)
  )

# PDFとして保存
ggsave("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/scatterplot_Ser2.pdf", width = 8, height = 6)
###################################################################
# 回帰直線の計算（Ser5_Aveを使用）
lm_red_Ser5_log2 <- lm(log2(Ser5_Ave) ~ Pol2_Ave, data = subset(df_highvar_no_df3, is_in_df3 == TRUE))
lm_black_Ser5_log2 <- lm(log2(Ser5_Ave) ~ Pol2_Ave, data = subset(df_highvar_no_df3, is_in_df3 == FALSE))

lm_red_Ser5 <- lm(Ser5_Ave ~ Pol2_Ave, data = subset(df_highvar_no_df3, is_in_df3 == TRUE))
lm_black_Ser5 <- lm(Ser5_Ave ~ Pol2_Ave, data = subset(df_highvar_no_df3, is_in_df3 == FALSE))

# 直線の式を文字列で作成
eq_red_Ser5_log2 <- paste("y = ", round(coef(lm_red_Ser5_log2)[2], 2), "x + ", round(coef(lm_red_Ser5_log2)[1], 2), sep = "")
eq_black_Ser5_log2 <- paste("y = ", round(coef(lm_black_Ser5_log2)[2], 2), "x + ", round(coef(lm_black_Ser5_log2)[1], 2), sep = "")

eq_red_Ser5 <- paste("y = ", round(coef(lm_red_Ser5)[2], 2), "x + ", round(coef(lm_red_Ser5)[1], 2), sep = "")
eq_black_Ser5 <- paste("y = ", round(coef(lm_black_Ser5)[2], 2), "x + ", round(coef(lm_black_Ser5)[1], 2), sep = "")

# x軸の範囲を指定
x_range <- range(df_highvar_no_df3$Pol2_Ave)

# log2(Ser5_Ave)を使用したプロット
ggplot(df_highvar_no_df3, aes(x = Pol2_Ave, y = log2(Ser5_Ave))) +
  geom_point(aes(color = is_in_df3), size = 2, alpha = 1) +  # 点の大きさを4に設定
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#2C3E50")) +  # TRUEは赤、FALSEは青
  geom_smooth(data = subset(df_highvar_no_df3, is_in_df3 == TRUE), 
              aes(x = Pol2_Ave, y = log2(Ser5_Ave)), 
              method = "lm", color = "red", se = FALSE) +  # 赤の回帰直線
  geom_smooth(data = subset(df_highvar_no_df3, is_in_df3 == FALSE), 
              aes(x = Pol2_Ave, y = log2(Ser5_Ave)), 
              method = "lm", color = "black", se = FALSE, 
              xlim = c(min(x_range), max(x_range) * 1.2)) +  # 黒の回帰直線を延長
  geom_text(aes(x = max(Pol2_Ave), y = log2(Ser5_Ave)[which.max(Pol2_Ave)], 
                label = eq_red_Ser5_log2), 
            color = "red", hjust = 1, vjust = 1) +  # 赤直線の式
  geom_text(aes(x = max(Pol2_Ave), y = log2(Ser5_Ave)[which.max(Pol2_Ave)], 
                label = eq_black_Ser5_log2), 
            color = "black", hjust = 1, vjust = 0) +  # 黒直線の式
  labs(
    x = "Pol2 Average",
    y = "log2(Ser5 Average)",  # 縦軸ラベルを変更
    title = "Pol2 vs log2(Ser5) Scatter Plot (High Variance Genes)",
    color = "In df3"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0.5)
  )

# log2(Ser5_Ave)でPDFとして保存
ggsave("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/scatterplot_Ser5_log2.pdf", width = 8, height = 6)

# Ser5_Aveを使用したプロット
ggplot(df_highvar_no_df3, aes(x = Pol2_Ave, y = Ser5_Ave)) +
  geom_point(aes(color = is_in_df3), size = 2, alpha = 1) +  # 点の大きさを4に設定
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#2C3E50")) +  # TRUEは赤、FALSEは青
  geom_smooth(data = subset(df_highvar_no_df3, is_in_df3 == TRUE), 
              aes(x = Pol2_Ave, y = Ser5_Ave), 
              method = "lm", color = "red", se = FALSE) +  # 赤の回帰直線
  geom_smooth(data = subset(df_highvar_no_df3, is_in_df3 == FALSE), 
              aes(x = Pol2_Ave, y = Ser5_Ave), 
              method = "lm", color = "black", se = FALSE, 
              xlim = c(min(x_range), max(x_range) * 1.2)) +  # 黒の回帰直線を延長
  geom_text(aes(x = max(Pol2_Ave), y = Ser5_Ave[which.max(Pol2_Ave)], 
                label = eq_red_Ser5), 
            color = "red", hjust = 1, vjust = 1) +  # 赤直線の式
  geom_text(aes(x = max(Pol2_Ave), y = Ser5_Ave[which.max(Pol2_Ave)], 
                label = eq_black_Ser5), 
            color = "black", hjust = 1, vjust = 0) +  # 黒直線の式
  labs(
    x = "Pol2 Average",
    y = "Ser5 Average",  # 縦軸ラベルを変更
    title = "Pol2 vs Ser5 Scatter Plot (High Variance Genes)",
    color = "In df3"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0.5)
  )

# Ser5_AveでPDFとして保存
ggsave("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/scatterplot_Ser5.pdf", width = 8, height = 6)
################################################################
# グラフ作成：Ser2_Aveの回帰直線なし、log2(Ser2_Ave)の回帰直線なし
ggplot(df_highvar_no_df3, aes(x = Pol2_Ave, y = Ser2_Ave)) +
  geom_point(aes(color = is_in_df3), size = 2, alpha = 1) +  # 点の大きさを4に設定
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#2C3E50")) +  # TRUEは赤、FALSEは青
  labs(
    x = "Pol2 Average",
    y = "Ser2 Average",  # 縦軸ラベルを変更
    title = "Pol2 vs Ser2 Scatter Plot (High Variance Genes) - No Regression Line",
    color = "In df3"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0.5)
  )

# Ser2_Aveの回帰直線なし、log2(Ser2_Ave)の回帰直線なしでPDF保存
ggsave("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/scatterplot_Ser2_no_regression.pdf", width = 8, height = 6)

# グラフ作成：log2(Ser2_Ave)の回帰直線なし
ggplot(df_highvar_no_df3, aes(x = Pol2_Ave, y = log2(Ser2_Ave))) +
  geom_point(aes(color = is_in_df3), size = 2, alpha = 1) +  # 点の大きさを4に設定
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#2C3E50")) +  # TRUEは赤、FALSEは青
  labs(
    x = "Pol2 Average",
    y = "log2(Ser2 Average)",  # 縦軸ラベルを変更
    title = "Pol2 vs log2(Ser2) Scatter Plot (High Variance Genes) - No Regression Line",
    color = "In df3"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0.5)
  )

# log2(Ser2_Ave)の回帰直線なしでPDF保存
ggsave("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/scatterplot_LogSer2_no_regression.pdf", width = 8, height = 6)


# グラフ作成：Ser5_Aveの回帰直線なし
ggplot(df_highvar_no_df3, aes(x = Pol2_Ave, y = Ser5_Ave)) +
  geom_point(aes(color = is_in_df3), size = 2, alpha = 1) +  # 点の大きさを4に設定
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#2C3E50")) +  # TRUEは赤、FALSEは青
  labs(
    x = "Pol2 Average",
    y = "Ser5 Average",  # 縦軸ラベルを変更
    title = "Pol2 vs Ser5 Scatter Plot (High Variance Genes) - No Regression Line",
    color = "In df3"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0.5)
  )

# Ser5_Aveの回帰直線なしでPDF保存
ggsave("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/scatterplot_Ser5_no_regression.pdf", width = 8, height = 6)


# グラフ作成：log2(Ser5_Ave)の回帰直線なし
ggplot(df_highvar_no_df3, aes(x = Pol2_Ave, y = log2(Ser5_Ave))) +
  geom_point(aes(color = is_in_df3), size = 2, alpha = 1) +  # 点の大きさを4に設定
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#2C3E50")) +  # TRUEは赤、FALSEは青
  labs(
    x = "Pol2 Average",
    y = "log2(Ser5 Average)",  # 縦軸ラベルを変更
    title = "Pol2 vs log2(Ser5) Scatter Plot (High Variance Genes) - No Regression Line",
    color = "In df3"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(hjust = 0.5)
  )

# log2(Ser5_Ave)の回帰直線なしでPDF保存
ggsave("C:/Documents/20241119Ser2Pol2scatterplot/bigwigAverageOverBed/scatterplot_LogSer5_no_regression.pdf", width = 8, height = 6)





