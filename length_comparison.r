# Rでの箱ひげ図を作成
setwd("C:/Documents/length_comparison")
# 必要なライブラリを読み込む
library(ggplot2)

# データの読み込み
included <- read.table("included_genes_subset.bed", header = FALSE)
excluded <- read.table("excluded_genes_subset.bed", header = FALSE)

# カラム名を設定
colnames(included) <- c("Length")
colnames(excluded) <- c("Length")
head(included)
head(excluded)
# それぞれのデータにラベルを追加
included$Group <- "Included"
excluded$Group <- "Excluded"

##################################################
# データの準備
included_mean <- mean(included$Length)
excluded_mean <- mean(excluded$Length)
included_median <- median(included$Length)
excluded_median <- median(excluded$Length)
u_test <- wilcox.test(included$Length, excluded$Length)

# 結果データフレーム（四捨五入＋有効数字3桁p値）
results <- data.frame(
  Metric = c("Mean", "Median"),
  Included = c(round(included_mean), round(included_median)),
  Excluded = c(round(excluded_mean), round(excluded_median)),
  p_value = c(NA, formatC(u_test$p.value, digits = 3, format = "e"))
)

# テーブル作成
library(gridExtra)
library(grid)

table_grob <- tableGrob(results, rows = NULL, theme = ttheme_minimal(base_size = 14))

# PDFに保存
pdf("C:/Documents/length_comparison/length_comparison_result.pdf", height = 5, width = 8)
grid.newpage()
grid.draw(table_grob)

# ヘッダー行（列名）とデータ1行目の間に横線を追加
# tableGrob では行0がヘッダー行に相当し、行1がデータの最初の行
y <- unit(1, "npc") - sum(table_grob$heights[1])
grid.lines(x = unit(c(0, 1), "npc"), y = y, gp = gpar(col = "black", lwd = 1.5))

dev.off()
##################################################
# 対数変換（底10）
included$LogLength <- log2(included$Length)
excluded$LogLength <- log2(excluded$Length)

# データを結合
combined_data <- rbind(included, excluded)

# 箱ひげ図を作成
plot <- ggplot(combined_data, aes(x = Group, y = LogLength)) +
  geom_boxplot(fill = c("lightblue", "lightgreen")) +
  theme_minimal() +
  labs(title = "Log-transformed Gene Length Comparison", y = "Log2(Gene Length)", x = "Group") +
  theme(text = element_text(size = 14))

# 画像を保存
ggsave("log_gene_length_comparison.pdf", plot = plot, height = 6, width = 3)

# ヒストグラムを作成（長さを横軸、遺伝子の個数の割合を縦軸）
plot <- ggplot(combined_data, aes(x = Length, fill = Group)) +
  geom_histogram(aes(y = ..density..), bins = 30, position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("lightblue", "lightgreen")) +
  theme_minimal() +
  labs(title = "Gene Length Distribution", x = "Gene Length", y = "Density (Proportion)") +
  theme(text = element_text(size = 14))

# 画像を保存
ggsave("gene_length_distribution_log2.pdf", plot = plot, height = 3, width = 10)

plott <- ggplot(combined_data, aes(x = log2(Length), fill = Group)) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.6) +
  scale_fill_manual(values = c("Included" = "lightgreen", "Excluded" = "lightblue")) +  # 大文字小文字を修正
  facet_wrap(~ Group, scales = "free_y") +
  theme_minimal() +
  labs(title = "Log2 Gene Length Distribution", x = "Log2(Gene Length)", y = "Density (Proportion)") +
  theme(text = element_text(size = 14), strip.text = element_text(size = 14))

ggsave("separated_gene_length_distribution_log2.pdf", plot = plott, height = 4, width = 10)

unique(combined_data$Group)

plott <- ggplot(combined_data, aes(x = log2(Length), fill = Group, color = Group)) +
  geom_density(alpha = 0.4, adjust = 1.5) +  # KDE曲線（adjustで滑らかさ調整）
  scale_fill_manual(values = c("Included" = "lightgreen", "Excluded" = "lightblue")) +
  scale_color_manual(values = c("Included" = "green", "Excluded" = "blue")) + 
  facet_wrap(~ Group, scales = "free_y") +
  theme_minimal() +
  labs(title = "Log2 Gene Length Density Distribution", x = "Log2(Gene Length)", y = "Density") +
  theme(text = element_text(size = 14), strip.text = element_text(size = 14))

ggsave("gene_length_density_distribution_log2.pdf", plot = plott, height = 4, width = 10)

# バイオリンプロットを使って分布を比較
p <- ggplot(combined_data, aes(x = Group, y = log2(Length), fill = Group)) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  scale_fill_manual(values = c("Included" = "lightblue", "Excluded" = "lightgreen")) +
  theme_minimal() +
  labs(title = "Log2-transformed Gene Length Violin Plot", x = "Group", y = "Log2(Gene Length)") +
  theme(text = element_text(size = 14))
ggsave("gene_length_violin_log2.pdf", plot = p, height = 10, width = 5)

# スワームプロット（データ点を個別に表示）
q <- ggplot(combined_data, aes(x = Group, y = log2(Length), color = Group)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.5) +
  scale_color_manual(values = c("Included" = "blue", "Excluded" = "green")) +
  theme_minimal() +
  labs(title = "Log2-transformed Gene Length Swarm Plot", x = "Group", y = "Log2(Gene Length)") +
  theme(text = element_text(size = 14))
ggsave("gene_length_swarm_log2.pdf", plot = q, height = 10, width = 5)
