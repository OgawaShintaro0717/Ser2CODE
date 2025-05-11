# プロット作成
g1 <- ggplot(y, aes(x = Distance_.microns.)) +
  geom_line(aes(y = AVIDIN, color = "AVIDIN"), size = 1.5) +
  geom_line(aes(y = COILIN, color = "COILIN"), size = 1.5)

# プロットの内部データを取得
plot_data <- ggplot_build(g1)

# 折れ線ごとの色情報を抽出
colors <- unique(plot_data$data[[2]]$colour)  # 1つ目のレイヤー (geom_line) の色

# 結果を確認
print(colors)
# プロットの内部データを取得
plot_data <- ggplot_build(g1)

# 折れ線ごとの色情報を抽出
line_colors <- unique(plot_data$data[[1]]$colour)

# 背景色を取得
background_color <- g1$theme$panel.background$fill
plot_background_color <- g1$theme$plot.background$fill

# 結果を表示
print("折れ線の色:")
print(line_colors)

print("プロットエリアの背景色:")
print(background_color)

print("全体の背景色:")
print(plot_background_color)
