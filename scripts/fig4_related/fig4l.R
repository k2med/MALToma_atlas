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
seurat_dc <- subset(
  seurat_myeloid,
  fine_cell_type %in% c("cDC_XCR1", "cDC_CLEC10A", "pDC_CLEC4C")
)
seurat_dc$fine_cell_type <- droplevels(seurat_dc$fine_cell_type)

# ---- build gene sets for AUCell scoring ----
manual_gene_sets <- list(
  Activated_DC = c('FSCN1',	'BIRC3',	'LAMP3',	'CCL19',	'LAD1',	'MARCKS',	'TNFAIP2',
                   'CCR7',	'CCL22',	'MARCKSL1',	'EBI3',	'TNFRSF11B',	'NUB1',	'INSM1',
                   'RAB9A',	'LY75',	'SIAH2',	'POGLUT1',	'KDM2B',	'MGLL',	'TXN',	'MLLT6',
                   'KIF2A',	'GRSF1',	'FAM49A',	'PLEKHG1',	'SOCS2',	'RFTN1',	'AC009812.4',
                   'BMP2K',	'NAV1',	'IL7R',	'ID2',	'CCL17',	'PPP1R9B',	'NRP2',	'TUBB6',
                   'ARNTL2',	'UVRAG',	'TXNDC11',	'MREG',	'BTG1',	'NDE1',	'SPG11',	'IL32',
                   'ERICH1',	'TBC1D4',	'NFKB1',	'GCSAM',	'BZW1'),
  Migratory_DC = c('GAL3ST',	'NUDT17',	'ITGB8',	'ADCY6',	'ENO2',	'IL15RA',	'SOCS2',
                   'IL15',	'STAP2',	'PHF24',	'ANKRD33B',	'INSM1',	'ANXA3',	'ARHGAP28',
                   'RNF115',	'ADORA2A',	'EXTL1',	'SPSB',	'SLC22A23',	'RABGAP1',	'GYG1',
                   'DAP',	'OGFR',	'GYG2',	'CCSER2',	'TMEM123',	'NET1',	'GPR52',	'SLCO5A1',
                   'FAH',	'CLU',	'PCGF5',	'SAMSN1',	'CDKN2B',	'BMP2K',	'ZC2HC1A',	'SERINC5',
                   'HIVEP1',	'CNR1',	'CNR2')
)

path_ids <- c(
  "Antigen processing and presentation [GOBP]",
  "Antigen processing and presentation of peptide antigen via mhc class i [GOBP]",
  "Antigen processing and presentation of peptide or polysaccharide antigen via mhc class ii [GOBP]"
)

pathway_gene_sets <- setNames(
  lapply(path_ids, function(pid) {
    unique(gene_sets[gene_sets$pathID_rename == pid, "geneID"])
  }),
  path_ids
)

gene_sets_list <- c(manual_gene_sets, pathway_gene_sets)

# ---- AUCell scoring on DCs ----
expr_counts   <- as.matrix(seurat_dc@assays$RNA@counts)
cell_rankings <- AUCell_buildRankings(expr_counts, nCores = 12, plotStats = TRUE)
cells_auc     <- AUCell_calcAUC(gene_sets_list, cell_rankings)

auc_score_mat <- AUCell::getAUC(cells_auc)  # gene sets x cells
stopifnot(identical(rownames(seurat_dc@meta.data), colnames(auc_score_mat)))
seurat_dc@meta.data <- cbind(seurat_dc@meta.data, t(auc_score_mat))

# ---- summarize by fine cell type and normalize ----
cell_type_score <- seurat_dc@meta.data[, c("fine_cell_type", names(gene_sets_list))] %>%
  dplyr::group_by(fine_cell_type) %>%
  dplyr::summarise(dplyr::across(.cols = where(is.numeric), .fns = median, na.rm = TRUE)) %>%
  tibble::column_to_rownames(var = "fine_cell_type")

cell_type_score_norm <- t(scale(cell_type_score))

# row splits should match the number/order of gene sets
row_split_vector <- c(1,1,2,2,2)

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

pdf(file.path(fig_dir, "fig4l.pdf"), width = 8, height = 4)
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
  height = unit(nrow(cell_type_score_norm) * 0.24 + 1 * 0.06, "cm"),
  row_split       = row_split_vector,
  top_annotation  = col_ann,
  heatmap_legend_param = heatmap_legend_param
)
dev.off()