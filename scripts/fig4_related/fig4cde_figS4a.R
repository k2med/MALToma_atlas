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
seurat_t <- readRDS(file.path(data_dir, "snrna_t_subset_seurat.rds"))

# ---- fine marker dotplot ----
Idents(seurat_t) <- seurat_t$fine_cell_type
t_fine_markers <- c(
  "CD8A", "CD4", "CCL5", "GZMK", "GZMA", "KLRG1",
  "GNLY", "GZMB", "KLRD1",
  "HAVCR2", "LAG3", "PDCD1",
  "FOXP3", "IL2RA", "TNFRSF4",
  "CXCL13", "TOX", "CXCR5",
  "TCF7", "CCR7", "SELL"
)

dp <- DotPlot(seurat_t, features = rev(t_fine_markers))
dotplot_df <- dp$data[, c("id", "features.plot", "pct.exp", "avg.exp.scaled")]
dotplot_df$id <- factor(dotplot_df$id, levels = levels(dotplot_df$id))

t_fine_dotplot <- ggplot(dotplot_df, aes(x = id, y = features.plot)) +
  geom_point(aes(color = avg.exp.scaled, size = pct.exp), shape = 19, stroke = 0.01) +
  labs(x = NULL, y = NULL) +
  scale_color_gradientn(colors = scale_col2) +
  scale_size(range = c(0, 2.4), limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  theme_bw(base_size = 5) +
  theme(
    text = element_text(size = 5),
    panel.grid.major = element_line(colour = "grey90", size = 0.116),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", size = 0.116),
    axis.line = element_line(colour = "black", size = 0.116),
    axis.ticks = element_line(colour = "black", size = 0.116),
    axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust = 1, size = 5),
    axis.text.y = element_text(colour = "black", size = 5, face = "italic"),
    legend.position = "right",
    legend.title = element_text(size = 5)
  ) +
  guides(
    size = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(stroke = 0.4)),
    fill = guide_colourbar(title.position = "top", title.hjust = 0.5)
  )

ggsave(
  filename = file.path(fig_dir, "fig4c.pdf"),
  plot = egg::set_panel_size(
    t_fine_dotplot,
    height = unit(length(t_fine_markers) / 6, "cm"),
    width  = unit(length(levels(dotplot_df$id)) / 6, "cm")
  ),
  width = 10, height = 10, units = "cm"
)

# ---- build Hallmark gene sets for AUCell scoring ----
gene_set_oxphos     <- gene_sets %>% dplyr::filter(pathID == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>% dplyr::pull(geneID)
gene_set_glycolysis <- gene_sets %>% dplyr::filter(pathID == "HALLMARK_GLYCOLYSIS") %>% dplyr::pull(geneID)
gene_set_fatty      <- gene_sets %>% dplyr::filter(pathID == "HALLMARK_FATTY_ACID_METABOLISM") %>% dplyr::pull(geneID)
gene_set_apoptosis  <- gene_sets %>% dplyr::filter(pathID == "HALLMARK_APOPTOSIS") %>% dplyr::pull(geneID)
gene_set_ifna       <- gene_sets %>% dplyr::filter(pathID == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>% dplyr::pull(geneID)
gene_set_ifng       <- gene_sets %>% dplyr::filter(pathID == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>% dplyr::pull(geneID)

gene_sets_list <- list(
  Naive          = c("CCR7", "IL7R", "LEF1", "SELL", "TCF7"),
  Cytotoxicity   = c("CST7", "CTSA", "GNLY", "GZMA", "GZMB", "IFNG", "KLRB1", "KLRD1", "KLRG1", "NKG7", "PRF1"),
  Exhaustion     = c("CTLA4", "ENTPD1", "HAVCR2", "LAG3", "PDCD1", "TIGIT", "TOX"),
  Treg           = c("CCR4", "CTLA4", "FOXP3", "ICOS", "IKZF2", "IKZF4", "IL10RA", "IL2RA", "TGFB1"),
  Proliferation  = c("AURKA", "BUB1", "CCNB1", "CCND1", "CCNE1", "DEK", "E2F1", "FEN1", "FOXM1",
                     "H2AFZ", "HMGB2", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MKI67", "MYBL2", "PCNA",
                     "PLK1", "TOP2A", "TYMS", "ZWINT"),
  Cytokine       = c("CSF1", "CSF2", "IL10RA", "IL16", "IL17RA", "IL18RAP", "IL6R", "IL21R", "IL2RA",
                     "IL2RB", "IL12RB2", "IL2RG", "IL32", "IL1R1", "IL1R2", "IL21", "IL26", "IL9R",
                     "ADAM10", "ADAM8", "ADAM12", "ADAM19", "METRNL", "CD70", "CXCL8", "XCL1", "XCL2",
                     "TGFB1", "TGFBR2", "TGFBR3"),
  Chemokine      = c("CCR4", "CCR5", "CCR6", "CCR7", "CCR8", "CXCR3", "CXCR4", "CXCR5", "CXCR6",
                     "CCL3", "CCL4", "CCL5", "CCL20", "CXCL13", "CXCL8", "XCL1", "XCL2"),
  Nfkb           = c("NFKB1", "NFKB2", "NFKBIA", "NFKBIB", "NFKBIZ", "CHUK", "IKBKB", "IKBKG", "REL", "RELA", "RELB"),
  IFNA           = gene_set_ifna,
  IFNG           = gene_set_ifng,
  Apoptosis      = gene_set_apoptosis,
  OXPHOS         = gene_set_oxphos,
  Glycolysis     = gene_set_glycolysis,
  Fatty_acid_metabolism = gene_set_fatty
)

# ---- AUCell scoring on T cells ----
expr_counts   <- as.matrix(seurat_t@assays$RNA@counts)
cell_rankings <- AUCell_buildRankings(expr_counts, nCores = 12, plotStats = TRUE)
cells_auc     <- AUCell_calcAUC(gene_sets_list, cell_rankings)

auc_score_mat <- AUCell::getAUC(cells_auc)  # gene sets x cells
stopifnot(identical(rownames(seurat_t@meta.data), colnames(auc_score_mat)))
seurat_t@meta.data <- cbind(seurat_t@meta.data, t(auc_score_mat))

# ---- summarize by fine cell type and normalize ----
cell_type_score <- seurat_t@meta.data[, c("fine_cell_type", names(gene_sets_list))] %>%
  dplyr::group_by(fine_cell_type) %>%
  dplyr::summarise(dplyr::across(.cols = where(is.numeric), .fns = median, na.rm = TRUE)) %>%
  tibble::column_to_rownames(var = "fine_cell_type")

cell_type_score_norm <- t(scale(cell_type_score))

row_split_vector <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3)
col_split_vector <- c(1, 1, 1, 2, 2, 2)

col_ann <- HeatmapAnnotation(
  celltype = colnames(cell_type_score_norm),
  col = list(celltype = fine_cell_type_col[colnames(cell_type_score_norm)]),
  simple_anno_size = unit(0.1, "cm"),
  annotation_name_gp = gpar(fontsize = 5)
)
heatmap_legend_param <- list(
  title_gp = gpar(fontsize = 5),
  labels_gp = gpar(fontsize = 5)
)

pdf(file.path(fig_dir, "fig4d.pdf"), width = 8, height = 4)
Heatmap(
  cell_type_score_norm,
  col = heatmap_col2,
  name = "Average score (Z-score)",
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  show_column_names = FALSE,
  rect_gp = gpar(col = "white", lwd = 0.33),
  show_heatmap_legend = TRUE,
  column_title_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 5),
  row_title_gp = gpar(fontsize = 5),
  row_names_gp = gpar(fontsize = 5),
  width  = unit(ncol(cell_type_score_norm) * 0.2 + 1 * 0.05, "cm"),
  height = unit(nrow(cell_type_score_norm) * 0.2 + 2 * 0.05, "cm"),
  row_split    = row_split_vector,
  column_split = col_split_vector,
  top_annotation = col_ann,
  heatmap_legend_param = heatmap_legend_param
)
dev.off()

# ---- CD8T subset boxplots (Cytotoxicity/Exhaustion) ----
meta_df_cd8t <- seurat_t@meta.data[seurat_t$fine_cell_type %in% c("CD8T_CCL5", "CD8T_GNLY", "CD8T_HAVCR2"), ]
meta_df_cd8t$fine_cell_type <- droplevels(meta_df_cd8t$fine_cell_type)

bx_input1 <- meta_df_cd8t %>%
  dplyr::select(subtype, fine_cell_type, Cytotoxicity, Exhaustion) %>%
  tidyr::pivot_longer(cols = c(Cytotoxicity, Exhaustion),
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

sig_boxplot1 <- wrap_plots(plots1, nrow = 3)
ggsave(filename = file.path(fig_dir, "fig4e.pdf"),
       plot = sig_boxplot1, width = 4, height = 6, units = "cm")

# ---- CD8T subset boxplots (Cytokine/Chemokine/IFN) ----
bx_input2 <- meta_df_cd8t %>%
  dplyr::select(subtype, fine_cell_type, Cytokine, Chemokine, IFNA, IFNG) %>%
  tidyr::pivot_longer(cols = c(Cytokine, Chemokine, IFNA, IFNG),
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

sig_boxplot2 <- wrap_plots(plots2, nrow = 3)
ggsave(filename = file.path(fig_dir, "figS4a.pdf"),
       plot = sig_boxplot2, width = 8, height = 6, units = "cm")