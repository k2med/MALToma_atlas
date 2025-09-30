# ---- packages ----
library(Seurat)
library(tidyverse)
library(AUCell)
library(patchwork)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))
gene_sets <- readRDS(file.path(src_dir, "msigdbr_gene_sets.rds"))

# ---- inputs ----
seurat_myeloid <- readRDS(file.path(data_dir, "snrna_myeloid_subset_seurat.rds"))
Idents(seurat_myeloid) <- seurat_myeloid$fine_cell_type

# ---- gene sets for AUCell ----
path_ids <- c(
  "Angiogenesis [HALLMARK]",
  "Neutrophil chemotaxis [GOBP]",
  "T cell chemotaxis [GOBP]",
  "Lymphocyte chemotaxis [GOBP]"
)

gene_sets_list <- setNames(
  lapply(path_ids, function(pid) unique(gene_sets[gene_sets$pathID_rename == pid, "geneID"])),
  path_ids
)

# ---- AUCell scoring ----
expr_counts   <- as.matrix(seurat_myeloid@assays$RNA@counts)
cell_rankings <- AUCell_buildRankings(expr_counts, nCores = 12, plotStats = TRUE)
cells_auc     <- AUCell_calcAUC(gene_sets_list, cell_rankings)

auc_score_mat <- AUCell::getAUC(cells_auc)  # gene sets x cells
stopifnot(identical(rownames(seurat_myeloid@meta.data), colnames(auc_score_mat)))
meta_df <- cbind(seurat_myeloid@meta.data, t(auc_score_mat))

# ---- radar coord ----
coord_radar <- function(theta = "x", start = 0, direction = 1) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto(
    "CoordRadar", CoordPolar,
    theta = theta, r = r, start = start,
    direction = sign(direction),
    is_linear = function(coord) TRUE
  )
}

# ---- constants for plotting ----
subtype_levels   <- c("LME_1", "LME_2", "LME_3")
original_labels  <- c("Mac_SPARC", "Mac_STAB1", "Mac_MARCO", "Mac_CHIT1",
                      "cDC_XCR1", "cDC_CLEC10A", "pDC_CLEC4C",
                      "Neu_G0S2", "Mast_CPA3")
short_labels     <- LETTERS[seq_along(original_labels)]
label_map        <- setNames(short_labels, original_labels)

make_safe_filename <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

# ---- loop over signatures ----
for (sig_id in names(gene_sets_list)) {
  radar_df <- meta_df[, c("subtype", "fine_cell_type", sig_id)]
  colnames(radar_df) <- c("subtype", "fine_cell_type", "value")
  
  radar_df$subtype <- factor(radar_df$subtype, levels = subtype_levels)
  
  radar_df <- radar_df %>%
    group_by(subtype, fine_cell_type) %>%
    summarise(mean = mean(value), sem = sd(value) / sqrt(n()), .groups = "drop")
  
  radar_df$fine_cell_type_short <- label_map[radar_df$fine_cell_type]
  radar_df$fine_cell_type_short <- factor(radar_df$fine_cell_type_short, levels = rev(short_labels))
  radar_df <- arrange(radar_df, fine_cell_type_short)
  
  radar_plot <- ggplot(
    radar_df,
    aes(x = fine_cell_type_short, y = mean, color = subtype, group = subtype, fill = subtype)
  ) +
    geom_point(size = 0.01) +
    geom_polygon(alpha = 0.3, linewidth = 0.116) +
    coord_radar() +
    theme_minimal() +
    scale_fill_manual(values = subtype_col1) +
    scale_color_manual(values = subtype_col2) +
    labs(x = NULL, y = NULL) +
    theme(
      text = element_text(size = 5),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.116),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(colour = "black", size = 5),
      axis.text.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 5),
      strip.background = element_blank()
    ) +
    scale_y_continuous(limits = c(0, max(radar_df$mean, na.rm = TRUE)))
  
  legend_df <- tibble(short = short_labels, label = names(label_map)) %>%
    mutate(short = factor(short, levels = rev(short_labels)))
  
  legend_plot <- ggplot(legend_df, aes(x = 1, y = short)) +
    geom_text(aes(label = paste(short, ":", label)), size = 5 / .pt, hjust = 0) +
    scale_y_discrete(limits = levels(legend_df$short)) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  merged_plot <- (radar_plot | legend_plot) + plot_layout(widths = c(1, 3))
  
  out_name <- sprintf("fig4m_%s.pdf", make_safe_filename(sig_id))
  ggsave(
    filename = file.path(fig_dir, out_name),
    plot = merged_plot,
    width = 10, height = 4, units = "cm"
  )
}