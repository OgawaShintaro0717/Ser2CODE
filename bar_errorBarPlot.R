#################################################################################
# 必要なライブラリを読み込む
library(ggplot2)
library(dplyr)

# データの作成
data <- read.csv("C:/Documents/20250102_7SK_RNAseq/matome.csv")

# CTと7SKのグループを分ける
data_CT <- data %>% filter(grepl("CT", sample))
data_7SK <- data %>% filter(grepl("7SK", sample))

# 各グループの平均値と標準誤差を計算
summary_CT <- data_CT %>%
  summarise(mean_CPM = mean(CPMreadCount), sd_CPM = sd(CPMreadCount), n = n())

summary_7SK <- data_7SK %>%
  summarise(mean_CPM = mean(CPMreadCount), sd_CPM = sd(CPMreadCount), n = n())

# 標準誤差（標準偏差 / sqrt(サンプル数)）を計算
summary_CT$se_CPM <- summary_CT$sd_CPM / sqrt(summary_CT$n)
summary_7SK$se_CPM <- summary_7SK$sd_CPM / sqrt(summary_7SK$n)

# 棒グラフのデータフレームを作成
plot_data <- data.frame(
  group = c("CT", "7SK"),
  mean = c(summary_CT$mean_CPM, summary_7SK$mean_CPM),
  se = c(summary_CT$se_CPM, summary_7SK$se_CPM)
)

# group を factor 型にして順序を指定
plot_data$group <- factor(plot_data$group, levels = c("CT", "7SK"))

# ggplotでエラー棒付きの棒グラフを描く
ggplot(plot_data, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.4) +  # 棒グラフを細くする
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +  # エラー棒
  labs(y = "Average CPM Read Count", x = "") +
  theme_minimal() +
  theme(legend.position = "none",  # 凡例を消す
        axis.title.x = element_blank(), # x軸タイトルを非表示
        axis.text.x = element_text(size = 12, face = "bold"),  # x軸のテキストを調整
        axis.title.y = element_text(size = 16, face = "bold")) +  # y軸タイトルを大きく
  scale_fill_manual(values = c("skyblue", "dodgerblue")) +  # 青系統の色を指定
  scale_x_discrete(labels = c("CT" = "siCT", "7SK" = "si7SK"))  # x軸ラベルを変更
