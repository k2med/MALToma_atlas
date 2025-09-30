# ---- packages ----
library(Seurat)
library(tidyverse)
library(ggplot2)
library(AUCell)
library(ComplexHeatmap)
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
seurat_fibro <- subset(
  seurat_stromal,
  fine_cell_type %in% c("Fibro_PI16", "Fibro_CCL19", "Fibro_LRRC15")
)
seurat_fibro$fine_cell_type <- droplevels(seurat_fibro$fine_cell_type)

# ---- build gene sets for AUCell scoring ----
path_ids <- c(
  "Antigen processing and presentation [GOBP]",
  "Lymphocyte chemotaxis [GOBP]",
  "Tnfa signaling via nfkb [HALLMARK]",
  "Extracellular matrix assembly [GOBP]",
  "Extracellular matrix constituent secretion [GOBP]",
  "Extracellular matrix organization [REACTOME]"
)

gene_sets_list <- setNames(
  lapply(path_ids, function(pid) {
    unique(gene_sets[gene_sets$pathID_rename == pid, "geneID"])
  }),
  path_ids
)

# ---- AUCell scoring on fibroblasts ----
expr_counts   <- as.matrix(seurat_fibro@assays$RNA@counts)
cell_rankings <- AUCell_buildRankings(expr_counts, nCores = 12, plotStats = TRUE)
cells_auc     <- AUCell_calcAUC(gene_sets_list, cell_rankings)

auc_score_mat <- AUCell::getAUC(cells_auc)  # gene sets x cells
stopifnot(identical(rownames(seurat_fibro@meta.data), colnames(auc_score_mat)))
seurat_fibro@meta.data <- cbind(seurat_fibro@meta.data, t(auc_score_mat))

# ---- summarize by fine cell type and normalize ----
cell_type_score <- seurat_fibro@meta.data[, c("fine_cell_type", names(gene_sets_list))] %>%
  dplyr::group_by(fine_cell_type) %>%
  dplyr::summarise(dplyr::across(.cols = where(is.numeric), .fns = median, na.rm = TRUE)) %>%
  tibble::column_to_rownames(var = "fine_cell_type")

cell_type_score_norm <- t(scale(cell_type_score))

col_ann <- HeatmapAnnotation(
  celltype = colnames(cell_type_score_norm),
  col = list(celltype = fine_cell_type_col[colnames(cell_type_score_norm)]),
  simple_anno_size   = unit(0.1, "cm"),
  annotation_name_gp = gpar(fontsize = 5)
)
heatmap_legend_param <- list(
  title_gp  = gpar(fontsize = 5),
  labels_gp = gpar(fontsize = 5)
)

pdf(file.path(fig_dir, "fig5b.pdf"), width = 8, height = 4)
Heatmap(
  cell_type_score_norm,
  col = heatmap_col2,
  name = "Average score (Z-score)",
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
  width  = unit(ncol(cell_type_score_norm) * 0.24, "cm"),
  height = unit(nrow(cell_type_score_norm) * 0.24, "cm"),
  top_annotation  = col_ann,
  heatmap_legend_param = heatmap_legend_param
)
dev.off()

# ---- Fibro subset boxplots (TNFa/ECM organization) ----
meta_df_fibro <- seurat_fibro@meta.data[seurat_fibro$fine_cell_type %in% c("Fibro_CCL19", "Fibro_LRRC15"), ]
meta_df_fibro$fine_cell_type <- droplevels(meta_df_fibro$fine_cell_type)

gene_sets_name1 <- c("Tnfa signaling via nfkb [HALLMARK]",
                     "Extracellular matrix organization [REACTOME]")

bx_input1 <- meta_df_fibro %>%
  dplyr::select(subtype, fine_cell_type, gene_sets_name1) %>%
  tidyr::pivot_longer(cols = gene_sets_name1,
                      names_to = "signature", values_to = "score")

plots1 <- bx_input1 %>%
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

sig_boxplot1 <- wrap_plots(plots1, nrow = 2)
ggsave(filename = file.path(fig_dir, "fig5g.pdf"),
       plot = sig_boxplot1, width = 4, height = 4, units = "cm")

# ---- Fibro subset boxplots (other gene sets) ----
gene_sets_name2 <- c("Antigen processing and presentation [GOBP]",
                     "Lymphocyte chemotaxis [GOBP]",
                     "Extracellular matrix assembly [GOBP]",
                     "Extracellular matrix constituent secretion [GOBP]")

bx_input2 <- meta_df_fibro %>%
  dplyr::select(subtype, fine_cell_type, gene_sets_name2) %>%
  tidyr::pivot_longer(cols = gene_sets_name2,
                      names_to = "signature", values_to = "score")

plots2 <- bx_input2 %>%
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

sig_boxplot2 <- wrap_plots(plots2, nrow = 2)
ggsave(filename = file.path(fig_dir, "figS4e.pdf"),
       plot = sig_boxplot2, width = 8, height = 4, units = "cm")
