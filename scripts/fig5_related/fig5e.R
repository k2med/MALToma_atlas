# ---- packages ----
library(Seurat)
library(tidyverse)
library(progeny)
library(ComplexHeatmap)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- inputs ----
seurat_stromal <- readRDS(file.path(data_dir, "snrna_stromal_subset_seurat.rds"))

# ---- subsets ----
seurat_fibro <- subset(
  seurat_stromal,
  fine_cell_type %in% c("Fibro_PI16", "Fibro_CCL19", "Fibro_LRRC15")
)
seurat_fibro$fine_cell_type <- droplevels(seurat_fibro$fine_cell_type)

# ---- PROGENy activity (as assay) ----
# Adds a "progeny" assay, then scale it
seurat_fibro <- progeny(
  seurat_fibro, scale = FALSE, organism = "Human", top = 500, perm = 1, return_assay = TRUE
) %>%
  ScaleData(assay = "progeny")

# ---- long-format activity with cell type ----
progeny_scores_df <- GetAssayData(seurat_fibro, slot = "scale.data", assay = "progeny") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  tidyr::pivot_longer(-cell, names_to = "pathway", values_to = "activity")

Idents(seurat_fibro) <- seurat_fibro$fine_cell_type
cells_clusters_df <- data.frame(
  cell = names(Idents(seurat_fibro)),
  cell_type = as.character(Idents(seurat_fibro)),
  stringsAsFactors = FALSE
)

progeny_scores_df <- progeny_scores_df %>%
  dplyr::inner_join(cells_clusters_df, by = "cell")

# ---- summarize mean activity per cell type ----
summarized_df <- progeny_scores_df %>%
  dplyr::group_by(pathway, cell_type) %>%
  dplyr::summarise(avg = mean(activity), .groups = "drop")

# ---- matrix (rows: pathway, cols: cell type) & z-score by pathway ----
wide_df <- summarized_df %>%
  tidyr::pivot_wider(names_from = pathway, values_from = avg) %>%
  tibble::column_to_rownames("cell_type")

# scale across cell types per pathway (row-wise after transpose)
scaled_mat <- t(scale(as.matrix(wide_df)))

# ---- subset pathways & cell types (order fixed) ----
pathways_keep <- c("JAK-STAT", "TNFa", "NFkB", "TGFb", "PI3K", "WNT")
celltypes_keep <- c("Fibro_PI16", "Fibro_CCL19", "Fibro_LRRC15")

heatmap_input <- scaled_mat[pathways_keep, celltypes_keep, drop = FALSE]

# ---- annotations & legend ----
col_ann <- HeatmapAnnotation(
  celltype = colnames(heatmap_input),
  col = list(celltype = fine_cell_type_col[colnames(heatmap_input)]),
  simple_anno_size   = unit(0.1, "cm"),
  annotation_name_gp = gpar(fontsize = 5)
)
heatmap_legend_param <- list(
  title_gp  = gpar(fontsize = 5),
  labels_gp = gpar(fontsize = 5)
)

# ---- heatmap ----
pdf(file.path(fig_dir, "fig5e.pdf"), width = 8, height = 4)
Heatmap(
  heatmap_input,
  col                  = heatmap_col2,
  name                 = "Average score (z-score)",
  cluster_columns      = FALSE,
  cluster_rows         = FALSE,
  show_row_dend        = FALSE,
  show_column_names    = FALSE,
  rect_gp              = gpar(col = "white", lwd = 0.33),
  show_heatmap_legend  = TRUE,
  column_title_gp      = gpar(fontsize = 5),
  column_names_gp      = gpar(fontsize = 5),
  row_title_gp         = gpar(fontsize = 5),
  row_names_gp         = gpar(fontsize = 5),
  width  = unit(ncol(heatmap_input) * 0.24, "cm"),
  height = unit(nrow(heatmap_input) * 0.24, "cm"),
  top_annotation       = col_ann,
  heatmap_legend_param = heatmap_legend_param
)
dev.off()