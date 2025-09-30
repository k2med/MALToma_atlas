# ---- packages ----
library(Seurat)
library(tidyverse)
library(AUCell)
library(ComplexHeatmap)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))
gene_sets <- readRDS(file.path(src_dir, "msigdbr_gene_sets.rds"))

# ---- inputs: main objects ----
seurat_myeloid <- readRDS(file.path(data_dir, "snrna_myeloid_subset_seurat.rds"))

# ---- subsets ----
seurat_mac <- subset(
  seurat_myeloid,
  fine_cell_type %in% c("Mac_SPARC", "Mac_STAB1", "Mac_MARCO", "Mac_CHIT1")
)
seurat_mac$fine_cell_type <- droplevels(seurat_mac$fine_cell_type)

# ---- build gene sets for AUCell scoring ----
manual_gene_sets <- list(
  M1_phenotype = c("IL23", "TNF", "CXCL9", "CXCL10", "CXCL11", "CD86", "IL1A", "IL1B", "IL6",
                   "CCL5", "IRF5", "IRF1", "CD40", "IDO1", "KYNU", "CCR7"),
  M2_phenotype = c("IL4R", "CCL4", "CCL13", "CCL20", "CCL17", "CCL18", "CCL22", "CCL24",
                   "LYVE1", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "EGF", "CTSA", "CTSB", "CTSC", "CTSD",
                   "TGFB1", "TGFB2", "TGFB3", "MMP14", "MMP19", "MMP9", "CLEC7A", "WNT7B", "FASL",
                   "TNFSF12", "TNFSF8", "CD276", "VTCN1", "MSR1", "FN1", "IRF4")
)

path_ids <- c(
  "Angiogenesis [HALLMARK]",
  "Neutrophil chemotaxis [GOBP]",
  "Phagocytosis [GOBP]",
  "Toll like receptor signaling pathway [KEGG]",
  "Antigen processing and presentation [GOBP]",
  "Regulation of t helper 1 cell differentiation [GOBP]",
  "Regulation of t helper 2 cell differentiation [GOBP]"
)

pathway_gene_sets <- setNames(
  lapply(path_ids, function(pid) {
    unique(gene_sets[gene_sets$pathID_rename == pid, "geneID"])
  }),
  path_ids
)

gene_sets_list <- c(manual_gene_sets, pathway_gene_sets)

# ---- AUCell scoring on macrophages ----
expr_counts   <- as.matrix(seurat_mac@assays$RNA@counts)
cell_rankings <- AUCell_buildRankings(expr_counts, nCores = 12, plotStats = TRUE)
cells_auc     <- AUCell_calcAUC(gene_sets_list, cell_rankings)

auc_score_mat <- AUCell::getAUC(cells_auc)  # gene sets x cells
stopifnot(identical(rownames(seurat_mac@meta.data), colnames(auc_score_mat)))
seurat_mac@meta.data <- cbind(seurat_mac@meta.data, t(auc_score_mat))

# ---- summarize by fine cell type and normalize ----
cell_type_score <- seurat_mac@meta.data[, c("fine_cell_type", names(gene_sets_list))] %>%
  dplyr::group_by(fine_cell_type) %>%
  dplyr::summarise(dplyr::across(.cols = where(is.numeric), .fns = median, na.rm = TRUE)) %>%
  tibble::column_to_rownames(var = "fine_cell_type")

cell_type_score_norm <- t(scale(cell_type_score))

# row splits should match the number/order of gene sets
row_split_vector <- c(1, 1, 2, 2, 3, 3, 4, 4, 4)

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

pdf(file.path(fig_dir, "fig4i.pdf"), width = 8, height = 4)
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
  height = unit(nrow(cell_type_score_norm) * 0.24 + 3 * 0.06, "cm"),
  row_split       = row_split_vector,
  top_annotation  = col_ann,
  heatmap_legend_param = heatmap_legend_param
)
dev.off()