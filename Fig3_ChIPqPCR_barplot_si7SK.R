# 必要なパッケージの読み込み
library(dplyr)
library(ggplot2)
library(tidyr)
# CSVを読み込んだデータフレームを表示
data <- read.csv("C:/Documents/qPCR_t_test_pair/Fig3si7SK/si7SK_AFF4.csv")
# データの先頭部分を表示（デフォルトで6行表示）
head(data)

# 同じ Primer と Condition に対する Value を平均化する
data_summarized <- data %>%
  group_by(Primer, Condition) %>%
  summarize(Value = mean(Value), .groups = "drop")

# Primerごとにデータを整形
data_wide <- data %>%
  spread(key = Condition, value = Value)

# 片側t検定を行い、p値を計算
p_values <- sapply(1:nrow(data_wide), function(i) {
  # si7SK が siCT より増加しているか減少しているかを判定
  if (mean(data_wide$si7SK[i]) > mean(data_wide$siCT[i])) {
    t_test_result <- t.test(data_wide$siCT[i], data_wide$si7SK[i], alternative = "less", paired = TRUE)
  } else {
    t_test_result <- t.test(data_wide$siCT[i], data_wide$si7SK[i], alternative = "greater", paired = TRUE)
  }
  return(t_test_result$p.value)
})

# p値に基づく記号の割り当て
significance <- ifelse(p_values < 0.01, "***",
                       ifelse(p_values < 0.05, "**",
                              ifelse(p_values < 0.1, "*",
                                     ifelse(p_values < 0.2, sprintf("p=%.3f", p_values), "NS"))))

# 結果をデータフレームにまとめる
results <- data.frame(
  Primer = data_wide$Primer,
  Mean_siCT = rowMeans(data_wide$siCT, na.rm = TRUE),
  Mean_si7SK = rowMeans(data_wide$si7SK, na.rm = TRUE),
  p_value = p_values,
  significance = significance
)

# 棒グラフの作成
p<- ggplot(results, aes(x = Primer, y = Mean_si7SK - Mean_siCT, fill = significance)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = significance), vjust = -0.5) +
  labs(y = "Difference (si7SK - siCT)", x = "Primer") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# グラフを保存（指定されたパス）
ggsave("C:/Documents/qPCR_t_test_pair/Fig3si7SK/si7SK_vs_siCT_plot_1.pdf", plot = p, width = 8, height = 6)
