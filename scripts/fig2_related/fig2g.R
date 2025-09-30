# ---- packages ----
library(GSVA)
library(GSEABase)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(grid)

# ---- paths & sources ----
data_dir  <- "../../data/bulk_rnaseq"
src_dir   <- "../source"
fig_dir   <- "."

bulk_expr_path     <- file.path(data_dir, "bulk_expression_matrix.txt")
bulk_clinical_path <- file.path(data_dir, "bulk_clinical_metadata.txt")
gmt_path           <- file.path(src_dir, "In_house_LME_signatures.gmt")

source(file.path(src_dir, "custom_colors.R"))

# ---- clinical data (tumor only) ----
clinical_df <- read.table(
  bulk_clinical_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)
clinical_tumor_df <- clinical_df[clinical_df$Group == "Tumor", , drop = FALSE]

# ---- load expression & GSVA ----
bulk_expr_mat <- read.table(
  bulk_expr_path, header = TRUE, row.names = 1, check.names = FALSE
)
bulk_expr_mat <- bulk_expr_mat[, rownames(clinical_tumor_df), drop = FALSE]
bulk_expr_mat <- as.matrix(log2(bulk_expr_mat + 1))

gene_sets <- getGmt(gmt_path, sep = "\t", geneIdType = SymbolIdentifier())
gsva_score_mat <- gsva(bulk_expr_mat, gene_sets, method = "gsva", kcdf = "Gaussian")

# ---- subtype factor aligned to GSVA columns ----
subtype_factor <- factor(clinical_tumor_df[colnames(gsva_score_mat), "LME_subtype"])
ord_idx <- order(subtype_factor)

# ---- marker set order (rows) grouped by lineage ----
celltype_groups <- list(
  B      = c("B_Mal", "B_Nor"),
  T      = c("CD8T_CCL5","CD8T_GNLY","CD8T_HAVCR2","CD4T_FOXP3","CD4T_CXCL13","CD4T_IL7R"),
  Myeloid= c("Mac_SPARC","Mac_STAB1","Mac_MARCO","Mac_CHIT1","cDC_XCR1","cDC_CLEC10A","pDC_CLEC4C","Neu_G0S2","Mast_CPA3"),
  Endo   = c("Endo_ACKR1","Endo_CXCL12","Endo_CCL21","Peri_PDGFRB"),
  Fibro  = c("Fibro_PI16","Fibro_CCL19","Fibro_LRRC15","FDC_CR2")
)
celltype_order <- unlist(celltype_groups, use.names = FALSE)

# subset & order rows (cell types) and columns (samples by subtype)
heatmap_input <- as.matrix(gsva_score_mat[celltype_order, ord_idx, drop = FALSE])

# ---- row scaling (z-score per signature) ----
heatmap_input_scaled <- t(apply(heatmap_input, 1, scale))
colnames(heatmap_input_scaled) <- colnames(heatmap_input)

# ---- annotations ----
col_ann <- HeatmapAnnotation(
  subtype = subtype_factor[ord_idx],
  col = list(subtype = subtype_col1),
  simple_anno_size = grid::unit(0.15, "cm"),
  annotation_name_gp = grid::gpar(fontsize = 5)
)

avg_per_subtype <- t(apply(heatmap_input_scaled, 1, function(x) {
  tapply(x, subtype_factor[ord_idx], mean)
}))

row_ann <- HeatmapAnnotation(
  which = "row",
  Average = anno_lines(
    avg_per_subtype,
    size  = grid::unit(1, "mm"),
    width = grid::unit(1, "cm"),
    pch   = 16,
    border = FALSE,
    axis_param = list(
      side   = "bottom",
      at     = c(min(avg_per_subtype), max(avg_per_subtype)),
      labels = c("Low", "High"),
      labels_rot = 45,
      gp = grid::gpar(fontsize = 5, lwd = 0.33)
    ),
    gp      = grid::gpar(col = subtype_col1, fill = subtype_col1),
    add_points = TRUE,
    pt_gp   = grid::gpar(col = subtype_col2, fill = subtype_col2)
  ),
  annotation_name_gp = grid::gpar(fontsize = 5)
)

# compute row_split vector by groups (same length as rows)
row_split_vec <- rep(seq_along(celltype_groups), lengths(celltype_groups))

# ---- draw heatmap ----
pdf(file.path(fig_dir, "fig2g.pdf"), width = 6, height = 6)
Heatmap(
  heatmap_input_scaled,
  col = heatmap_col1,
  name = "GSVA score",
  cluster_columns = FALSE,
  cluster_rows    = FALSE,
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  width  = grid::unit(4, "cm"),
  height = grid::unit(5, "cm"),
  column_title_gp = grid::gpar(fontsize = 5),
  column_names_gp = grid::gpar(fontsize = 5),
  row_title_gp    = grid::gpar(fontsize = 5),
  row_names_gp    = grid::gpar(fontsize = 5),
  row_split       = row_split_vec,
  column_split    = subtype_factor[ord_idx],
  top_annotation  = col_ann,
  right_annotation = row_ann
)
dev.off()