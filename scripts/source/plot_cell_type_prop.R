plot_cell_type_proportion <- function(seurat_obj,
                                      all_seurat_meta,
                                      cell_type_col   = "fine_cell_type",
                                      group_col       = "subtype",
                                      orig_col        = "orig.ident",
                                      cell_type_colors,
                                      filename        = "celltype_proportion.pdf") {
  
  # bar data: group-level composition
  bar_input <- seurat_obj@meta.data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_col, cell_type_col)))) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) %>%
    dplyr::mutate(relative_freq = n / sum(n))
  
  bar_plot <- ggplot(
    data = bar_input,
    aes(x = .data[[group_col]],
        y = relative_freq,
        fill = .data[[cell_type_col]],
        stratum = .data[[cell_type_col]],
        alluvium = .data[[cell_type_col]])
  ) +
    geom_col(position = "fill", width = 0.5) +
    ggalluvial::geom_flow(width = 0.5, fill = "white", knot.pos = 0,
                          linewidth = 0.116, color = "black", linetype = "dashed") +
    scale_fill_manual(values = cell_type_colors) +
    theme(
      text = element_text(size = 5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.116),
      axis.ticks = element_line(colour = "black", linewidth = 0.116),
      axis.text.x = element_text(colour = "black", size = 5, angle = 60, vjust = 1, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 5),
      legend.position = "none"
    ) +
    labs(x = NULL, y = "Proportion (%)", title = NULL) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.5),
      labels = function(x) x * 100
    )
  
  # line data: per-sample composition
  line_input <- seurat_obj@meta.data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(orig_col, cell_type_col)))) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(orig_col))) %>%
    dplyr::mutate(relative_freq = n / sum(n)) %>%
    dplyr::ungroup() %>%
    tidyr::complete(tidyr::nesting(!!rlang::sym(orig_col)), !!rlang::sym(cell_type_col),
                    fill = list(n = 0, relative_freq = 0)) %>%
    dplyr::left_join(all_seurat_meta[, c(orig_col, group_col)], by = orig_col) %>%
    dplyr::distinct()
  
  # ensure group order matches stacked bars
  line_input[[group_col]] <- factor(line_input[[group_col]],
                                    levels = unique(bar_input[[group_col]]))
  
  line_plot <- ggpubr::ggline(
    line_input,
    x = group_col,
    y = "relative_freq",
    color = cell_type_col,
    size = 0.116,
    point.size = 0.01,
    shape = 16,
    add = "mean_se"
  )
  
  # add pairwise Wilcoxon annotations per cell type (stacked vertically)
  y_set <- 0.95
  for (ct in unique(seurat_obj@meta.data[[cell_type_col]])) {
    line_plot <- line_plot +
      ggpubr::stat_pwc(
        data = line_input[line_input[[cell_type_col]] == ct, ],
        aes(group = .data[[group_col]]),
        method = "wilcox.test",
        p.adjust.method = "BH",
        label = "p.adj.signif",
        vjust = 0.5,
        step.increase = 0.05,
        hide.ns = TRUE,
        label.size = 5 * 25.4 / 72,
        tip.length = 0,
        size = 0.116,
        color = cell_type_colors[ct],
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
      axis.text.x = element_text(colour = "black", size = 5, angle = 60, vjust = 1, hjust = 1),
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
  
  final_plot <- (bar_plot | line_plot) + plot_layout(widths = c(1, 1.5))
  
  ggsave(filename = file.path(fig_dir, filename),
         plot = final_plot, width = 5, height = 3.2, units = "cm")
}