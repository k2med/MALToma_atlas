# ---- packages ----
library(Seurat)
library(tidyverse)
library(ggplot2)
library(AUCell)
library(ggpubr)
library(patchwork)
library(egg)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))
gene_sets <- readRDS(file.path(src_dir, "msigdbr_gene_sets.rds"))

# ---- inputs: main objects ----
seurat_stromal <- readRDS(file.path(data_dir, "snrna_stromal_subset_seurat.rds"))

# ---- subsets ----
seurat_endo <- subset(
  seurat_stromal,
  fine_cell_type %in% c("Endo_ACKR1", "Endo_CXCL12", "Endo_CCL21")
)
seurat_vec <- subset(
  seurat_endo,
  fine_cell_type %in% c("Endo_ACKR1", "Endo_CXCL12")
)
seurat_endo$fine_cell_type <- droplevels(seurat_endo$fine_cell_type)
seurat_vec$fine_cell_type <- droplevels(seurat_vec$fine_cell_type)

# ---- gene sets for AUCell (vascular biology programs) ----
# These pathway names must match 'pathID_rename' in the msigdbr_gene_sets.rds
path_ids <- c(
  "Angiogenesis [HALLMARK]",
  "Hypoxia [HALLMARK]",
  "Cellular response to vascular endothelial growth factor stimulus [GOBP]",
  "Vascular endothelial cell proliferation [GOBP]",
  "Vascular endothelial growth factor signaling pathway [GOBP]",
  "Vasculature development [GOBP]"
)

gene_sets_list <- setNames(
  lapply(path_ids, function(pid) {
    unique(gene_sets[gene_sets$pathID_rename == pid, "geneID"])
  }),
  path_ids
)

# ---- AUCell scoring on vascular endothelial cells (Endo_ACKR1 / Endo_CXCL12) ----
expr_counts   <- as.matrix(seurat_vec@assays$RNA@counts)
cell_rankings <- AUCell_buildRankings(expr_counts, nCores = 12, plotStats = TRUE)
cells_auc     <- AUCell_calcAUC(gene_sets_list, cell_rankings)

auc_score_mat <- AUCell::getAUC(cells_auc)  # signatures x cells
stopifnot(identical(rownames(seurat_vec@meta.data), colnames(auc_score_mat)))
meta_df <- cbind(seurat_vec@meta.data, t(auc_score_mat))

# ---- helper: polar radar coordinate ----
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

# ---- constants for radar plots ----
# Short letter tags for pathways used only as compact labels on the circle
original_labels <- path_ids
short_labels    <- LETTERS[seq_along(original_labels)]
label_map       <- setNames(short_labels, original_labels)

# ---- Radar: by cell types ----
# Aggregate AUCell AUC by fine_cell_type across cells for each pathway
radar_df <- meta_df %>%
  dplyr::select(fine_cell_type, all_of(path_ids)) %>%
  tidyr::pivot_longer(cols = all_of(path_ids),
                      names_to = "pathway", values_to = "value") %>%
  dplyr::group_by(fine_cell_type, pathway) %>%
  dplyr::summarise(mean = mean(value), sem = sd(value) / sqrt(dplyr::n()), .groups = "drop")

# Map long pathway names to short letters for a clean circular axis
radar_df$pathway_short <- label_map[radar_df$pathway]
radar_df$pathway_short <- factor(radar_df$pathway_short, levels = rev(short_labels))
radar_df <- dplyr::arrange(radar_df, pathway_short)

# Radar plot: each cell type forms a polygon across pathway spokes
radar_plot <- ggplot(
  radar_df,
  aes(x = pathway_short, y = mean, color = fine_cell_type, group = fine_cell_type, fill = fine_cell_type)
) +
  geom_point(size = 0.01) +
  geom_polygon(alpha = 0.3, linewidth = 0.116) +
  coord_radar() +
  theme_minimal() +
  scale_fill_manual(values = fine_cell_type_col) +
  scale_color_manual(values = c("#bbba11", "#5bad29")) +  # customize if needed
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

# Add a legend panel mapping short letters back to full pathway names
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

ggsave(
  filename = file.path(fig_dir, "fig5j_top.pdf"),
  plot = merged_plot,
  width = 10, height = 4, units = "cm"
)

# ---- Radar: by LME subtype ----
# Same idea but aggregate by LME subtype across cells (in the same two VE subtypes)
radar_df <- meta_df %>%
  dplyr::select(subtype, all_of(path_ids)) %>%
  tidyr::pivot_longer(cols = all_of(path_ids),
                      names_to = "pathway", values_to = "value") %>%
  dplyr::group_by(subtype, pathway) %>%
  dplyr::summarise(mean = mean(value), sem = sd(value) / sqrt(dplyr::n()), .groups = "drop")

radar_df$pathway_short <- label_map[radar_df$pathway]
radar_df$pathway_short <- factor(radar_df$pathway_short, levels = rev(short_labels))
radar_df <- dplyr::arrange(radar_df, pathway_short)

radar_plot <- ggplot(
  radar_df,
  aes(x = pathway_short, y = mean, color = subtype, group = subtype, fill = subtype)
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

ggsave(
  filename = file.path(fig_dir, "fig5j_bottom.pdf"),
  plot = merged_plot,
  width = 10, height = 4, units = "cm"
)

# ---- LEC (Endo_CCL21) signatures ----
path_ids <- c("Lymphocyte chemotaxis [GOBP]", "Interferon signaling [REACTOME]")

gene_sets_list <- setNames(
  lapply(path_ids, function(pid) unique(gene_sets[gene_sets$pathID_rename == pid, "geneID"])),
  path_ids
)

# AUCell on all endothelial cells (to include Endo_CCL21 cells)
expr_counts   <- as.matrix(seurat_endo@assays$RNA@counts)
cell_rankings <- AUCell_buildRankings(expr_counts, nCores = 12, plotStats = TRUE)
cells_auc     <- AUCell_calcAUC(gene_sets_list, cell_rankings)

auc_score_mat <- AUCell::getAUC(cells_auc)
stopifnot(identical(rownames(seurat_endo@meta.data), colnames(auc_score_mat)))
seurat_endo@meta.data <- cbind(seurat_endo@meta.data, t(auc_score_mat))

# Filter to LEC cells, long-format for ggplot split
meta_df_lec <- seurat_endo@meta.data[seurat_endo$fine_cell_type == "Endo_CCL21", ]
meta_df_lec$fine_cell_type <- droplevels(meta_df_lec$fine_cell_type)

bx_input <- meta_df_lec %>%
  dplyr::select(subtype, fine_cell_type, all_of(path_ids)) %>%
  tidyr::pivot_longer(cols = all_of(path_ids),
                      names_to = "signature", values_to = "score")

plots <- bx_input %>%
  split(interaction(.$signature, .$fine_cell_type)) %>%
  lapply(function(df) {
    ggplot(df, aes(x = subtype, y = score, fill = subtype)) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.116, width = 0.6) +
      scale_fill_manual(values = subtype_col1) +
      theme_minimal(base_size = 5) +
      theme(
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.116),
        axis.ticks = element_line(colour = "black", size = 0.116),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 5),
        legend.position = "none",
        plot.title = element_text(size = 5, hjust = 0.5)
      ) +
      ggtitle(paste(unique(df$signature), unique(df$fine_cell_type))) +
      ggpubr::stat_pwc(
        aes(group = subtype),
        label.size = 5 * 25.4 / 72,
        tip.length = 0,
        size = 0.116,
        hide.ns = TRUE,
        vjust = 0.5,
        method = "wilcox.test",
        p.adjust.method = "BH",
        label = "p.adj.signif"
      ) +
      coord_cartesian(clip = "off") +
      labs(x = "", y = "")
  })

sig_boxplot <- wrap_plots(plots, nrow = 2)
ggsave(
  filename = file.path(fig_dir, "figS4g.pdf"),
  plot = sig_boxplot, width = 2.5, height = 4, units = "cm"
)