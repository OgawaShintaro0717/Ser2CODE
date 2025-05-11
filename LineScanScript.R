setwd("D:/molecularbiology_20230223/IFData/20240725_ogawa_Turbo-coilin")
library("tidyverse")
library("ggplot2")
library("ggsci")

head(y)
g1 <- ggplot(y, aes(x = Distance_.microns. , y = Intensity))
g1 <- g1 + geom_line(aes(y = CTD, color = "CTD"))
g1 <- g1 + geom_line(aes(y = COILIN, colour = "COILIN"))
plot(g1)

y <- read.csv("Values_Merge.csv")
head(y)
g1 <- ggplot(y, aes(x = Distance_.microns. , y = Intensity)) +
          geom_line(aes(y = AVIDIN, color = "AVIDIN" ),size = 1.5) +
          geom_line(aes(y = COILIN, color = "COILIN" ),size = 1.5)
plot(g1)


g1 <- ggplot(y, aes(x = Distance_.microns. , y = Intensity)) +
  geom_line(aes(y = COILIN, color = "COILIN" ),size = 1.5)+
  geom_line(aes(y = AVIDIN, color = "AVIDIN" ),size = 1.5) +

plot(g1)


# プロット作成
g1 <- ggplot(y, aes(x = Distance_.microns.)) +
  geom_line(aes(y = COILIN, color = "COILIN"), size = 1.5) +
  geom_line(aes(y = AVIDIN, color = "AVIDIN"), size = 1.5) +
  scale_color_manual(values = c("COILIN" = "#F8766D", "AVIDIN" = "#00BFC4")) +  # より落ち着いた赤と緑
  labs(x = "Distance (microns)", y = "Signal Intensity", color = "Signal Type") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "#F5F5F5", color = NA),  # 背景を薄い灰色に（白に近い）
    panel.grid.major = element_line(color = "#E0E0E0"),  # メジャーグリッド線を薄い灰色に
    panel.grid.minor = element_line(color = "#E0E0E0", linetype = "dotted"),  # マイナーグリッド線を薄い灰色の点線に
    legend.background = element_rect(fill = "white", color = NA),  # 凡例の背景を薄い灰色に
    legend.key = element_rect(fill = "#F5F5F5", color = NA)  # 凡例キーの背景も統一
  )

# プロット表示
plot(g1)


# PDFに保存（縦4 x 横7のサイズ）
ggsave("line_scan_plot.pdf", plot = g1, width = 7, height = 4, units = "in")
