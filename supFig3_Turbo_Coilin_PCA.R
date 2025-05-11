data <-read.csv("E:/molecularbiology/my_Article/supportFigure3/TurboID_coilin/TURBOCOILIN4PCA.csv")
head(data)
library(ggplot2)

# データ読み込み
data <- read.csv("E:/molecularbiology/my_Article/supportFigure3/TurboID_coilin/TURBOCOILIN4PCA.csv", row.names = 1)

# データの転置
data_t <- t(data)

# 分散ゼロの列を削除
constant_columns <- apply(data_t, 2, var) == 0
data_t_clean <- data_t[, !constant_columns]

# PCA解析
pca_result <- prcomp(data_t_clean, scale. = TRUE)
pca_data <- as.data.frame(pca_result$x)

# グループ情報を追加
pca_data$group <- rownames(pca_data)  # サンプル名をグループとして使用

# PC1 vs PC2のプロット（点の大きさを調整）
p2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 5) +  # 点の大きさを調整
  coord_fixed(ratio = 1.5) +  # 軸比を変更
  labs(title = "PCA Plot (PC1 vs PC2)", x = "PC1", y = "PC2") +
  theme_minimal()

# PDFとして保存
pdf_path <- "E:/molecularbiology/my_Article/supportFigure3/TurboID_coilin/PCA_TurboID_Coilin.pdf"
ggsave(filename = pdf_path, plot = p2, width = 6, height = 6, device = cairo_pdf)

# 確認用にプロットを表示
print(p2)

history(max.show=Inf)

history(max.show=Inf)

ls -lt ~/.Rproj.user
file.show(".Rhistory")


