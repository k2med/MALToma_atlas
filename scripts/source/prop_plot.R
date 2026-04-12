plot_cell_type_proportion <- function(seurat_obj, 
                                      all_seurat_meta, 
                                      cell_type_col = "fine_cell_type",
                                      group_col = "subtype",
                                      orig_col = "orig.ident",
                                      cell_type_colors,
                                      filename = "celltype_proportion.pdf") {
  require(dplyr)
  require(ggplot2)
  require(ggalluvial)
  require(ggpubr)
  require(tidyr)
  require(patchwork)
  
  # 柱状图数据
  bar_input <- seurat_obj@meta.data %>%
    group_by(across(all_of(c(group_col, cell_type_col)))) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(across(all_of(group_col))) %>%
    mutate(relative_freq = n / sum(n))
  
  # 柱状图
  bar_plot <- ggplot(
    data = bar_input,
    aes_string(x = group_col, y = "relative_freq", fill = cell_type_col,
               stratum = cell_type_col, alluvium = cell_type_col)) +
    geom_col(position = "fill", width = 0.5) +
    geom_flow(width = 0.5, fill = 'white', knot.pos = 0, size = 0.116,
              color = 'black', linetype = "dashed") +
    scale_fill_manual(values = cell_type_colors) +
    theme(
      text = element_text(size = 5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size = 0.116),
      axis.ticks = element_line(colour = "black", size = 0.116),
      axis.text.x = element_text(colour = "black", size = 5,
                                 angle = 60, vjust = 1, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 5),
      legend.position = "none"
    ) +
    labs(x = "", y = "Proportion (%)", title = '') +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.5),
      labels = function(x) x * 100
    )
  
  # 折线图数据
  line_input <- seurat_obj@meta.data %>%
    group_by(across(all_of(c(orig_col, cell_type_col)))) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(across(all_of(orig_col))) %>%
    mutate(relative_freq = n / sum(n)) %>%
    ungroup() %>%
    complete(nesting(!!sym(orig_col)), !!sym(cell_type_col),
             fill = list(n = 0, relative_freq = 0)) %>%
    left_join(all_seurat_meta[, c(orig_col, group_col)], by = orig_col) %>%
    distinct()
  
  # 设置 factor 顺序
  line_input[[group_col]] <- factor(line_input[[group_col]],
                                    levels = unique(bar_input[[group_col]]))
  
  # 折线图
  line_plot <- ggline(
    line_input,
    x = group_col,
    y = "relative_freq",
    color = cell_type_col,
    size = 0.116,
    point.size = 0.01,
    shape = 16,
    add = "mean_se"
  )
  
  y_set <- 0.95
  for (cell_type in unique(seurat_obj@meta.data[[cell_type_col]])) {
    line_plot <- line_plot +
      stat_pwc(
        data = line_input[line_input[[cell_type_col]] == cell_type, ],
        aes_string(group = group_col),
        method = "wilcox.test",
        p.adjust.method = "BH",
        label = "p.adj.signif",
        vjust = 0.5,
        step.increase = 0.05,
        hide.ns = TRUE,
        label.size = 5 * 25.4 / 72,
        tip.length = 0,
        size = 0.116,
        color = cell_type_colors[cell_type],
        y.position = y_set
      )
    y_set <- y_set - 0.1
  }
  
  line_plot <- line_plot +
    scale_color_manual(values = cell_type_colors) +
    theme(
      text = element_text(size = 5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size = 0.116),
      axis.ticks = element_line(colour = "black", size = 0.116),
      axis.text.x = element_text(colour = "black", size = 5,
                                 angle = 60, vjust = 1, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 5),
      legend.position = "none"
    ) +
    labs(x = "", y = "", title = '') +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.5),
      labels = function(x) x * 100
    ) +
    coord_cartesian(clip = "off")
  
  # 合并
  final_plot <- (bar_plot | line_plot) +
    plot_layout(widths = c(1, 1.5))
  
  ggsave(filename = filename, plot = final_plot, width = 5, height = 3.2, units = 'cm')
}