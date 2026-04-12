# ---- packages ----
library(GSVA)
library(GSEABase)
library(NMF)
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
bulk_expr_mat <- bulk_expr_mat[, rownames(clinical_tumor_df)]
bulk_expr_mat <- as.matrix(log2(bulk_expr_mat + 1))

gene_sets <- getGmt(gmt_path, sep = "\t", geneIdType = SymbolIdentifier())

gsva_score_mat <- gsva(bulk_expr_mat, gene_sets, method = "gsva", kcdf = "Gaussian")

# ---- normalization for NMF ----
min_max_normalize <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

nmf_input_mat <- t(apply(gsva_score_mat, 1, min_max_normalize))

# ---- NMF rank survey ----
nmf_estim <- nmf(nmf_input_mat, 2:7, nrun = 100, seed = 421, method = "snmf/r")

plot(nmf_estim)
ggsave(
  filename = file.path(fig_dir, "figS4b.pdf"),
  width = 18, height = 12, units = "cm"
)

# ---- final NMF at chosen rank ----
nmf_rank   <- 3
nmf_result <- nmf(
  nmf_input_mat,
  rank  = nmf_rank,
  nrun  = 100,
  seed  = 421,
  method = "snmf/r"
)

# Note:
# The naming order here was retained to remain consistent with the initial version for matching purposes.
# In the revised version, the order was rearranged following the reviewer’s comment.
# However, this difference in naming order does not affect the analysis workflow or the final results.

nmf_cluster <- NMF::predict(nmf_result)
lme_map <- c("1" = "2", "2" = "3", "3" = "1")

subtype_factor <- factor(
  paste0("LME_", lme_map[as.character(nmf_cluster)])
)

# ---- heatmap input (order by subtype) ----
ord_idx <- order(subtype_factor)
heatmap_input_mat <- gsva_score_mat[, ord_idx, drop = FALSE]
heatmap_scaled_mat <- t(apply(heatmap_input_mat, 1, scale))
colnames(heatmap_scaled_mat) <- colnames(heatmap_input_mat)

# ---- column annotations ----
# Subtype (ordered)
subtype_ordered <- subtype_factor[ord_idx]

clinical_tumor_df <- clinical_tumor_df[colnames(heatmap_scaled_mat), ]

age_bin_vec <- dplyr::case_when(
  clinical_tumor_df$Age >= 80 ~ ">=80",
  clinical_tumor_df$Age >= 70 ~ "70-79",
  clinical_tumor_df$Age >= 60 ~ "60-69",
  clinical_tumor_df$Age >= 50 ~ "50-59",
  clinical_tumor_df$Age >= 40 ~ "40-49",
  clinical_tumor_df$Age <  40 ~ "<40",
  TRUE ~ NA_character_
)

col_ann <- HeatmapAnnotation(
  Subtype   = subtype_ordered,
  Sex       = clinical_tumor_df$Sex,
  Age       = age_bin_vec,
  Stage     = clinical_tumor_df$Ann_Arbor_stage,
  Fusion    = clinical_tumor_df$BIRC3_MALT_fusion,
  Mutation  = clinical_tumor_df$TNFAIP3_mutation,
  col = list(
    Subtype   = subtype_col1,
    Sex       = sex_col,
    Age       = age_col,
    Stage     = stage_col,
    Fusion    = fusion_col,
    Mutation  = mutation_col
  ),
  simple_anno_size   = unit(0.15, "cm"),
  annotation_name_gp = gpar(fontsize = 5)
)

# ---- row-side trend per subtype (lines) ----
# average value per gene across each subtype split
avg_per_subtype_mat <- t(apply(heatmap_scaled_mat, 1, function(x) {
  tapply(x, subtype_ordered, mean, na.rm = TRUE)
}))

# ensure consistent column order
avg_per_subtype_mat <- avg_per_subtype_mat[, levels(subtype_ordered), drop = FALSE]

row_ann <- HeatmapAnnotation(
  which  = "row",
  Average = anno_lines(
    avg_per_subtype_mat,
    size       = unit(1, "mm"),
    width      = unit(1, "cm"),
    pch        = 16,
    border     = FALSE,
    axis_param = list(
      side   = "bottom",
      at     = c(min(avg_per_subtype_mat, na.rm = TRUE),
                 max(avg_per_subtype_mat, na.rm = TRUE)),
      labels = c("Low", "High"),
      labels_rot = 45,
      gp     = gpar(fontsize = 5, lwd = 0.33)
    ),
    gp          = gpar(col = subtype_col1, fill = subtype_col1),
    add_points  = TRUE,
    pt_gp       = gpar(col = subtype_col2, fill = subtype_col2)
  ),
  annotation_name_gp = gpar(fontsize = 5)
)

# ---- heatmap ----
pdf(file.path(fig_dir, "fig1g.pdf"), width = 6, height = 6)
Heatmap(
  heatmap_scaled_mat,
  col                 = heatmap_col1,
  name                = "GSVA score",
  cluster_columns     = FALSE,
  cluster_rows        = TRUE,
  row_dend_gp         = gpar(lwd = 0.33),
  row_dend_reorder    = FALSE,
  show_column_names   = FALSE,
  show_heatmap_legend = TRUE,
  width               = unit(4, "cm"),
  height              = unit(5, "cm"),
  column_title_gp     = gpar(fontsize = 5),
  column_names_gp     = gpar(fontsize = 5),
  row_title_gp        = gpar(fontsize = 5),
  row_names_gp        = gpar(fontsize = 5),
  column_split        = subtype_ordered,
  top_annotation      = col_ann,
  right_annotation    = row_ann
)
dev.off()
