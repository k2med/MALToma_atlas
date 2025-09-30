# ---- packages ----
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(egg)
library(plot3D)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- all cells ----
All_cell_seurat <- readRDS(file.path(data_dir, "snrna_all_cells_seurat.rds"))

Use_cell_seurat <- subset(
  All_cell_seurat,
  fine_cell_type %in% c("Fibro_TCF21", "AT1_AGER", "AT2_SFTPC", "Hepa_ALB", "Endo_HPGD"),
  invert = TRUE
)
Use_cell_seurat$fine_cell_type <- droplevels(Use_cell_seurat$fine_cell_type)

# ---- per-sample cell-type proportions ----
ft_by_sample <- table(Use_cell_seurat$fine_cell_type, Use_cell_seurat$orig.ident)
fine_cell_type_prop <- t(prop.table(ft_by_sample, margin = 2))
fine_cell_type_prop <- as.matrix(fine_cell_type_prop)

# ---- PCA ----
PCA_prcomp <- prcomp(fine_cell_type_prop, center = TRUE, scale. = TRUE)
PCA_pred   <- predict(PCA_prcomp)

# map each sample (orig.ident) to its subtype from metadata
sample_subtype <- vapply(
  rownames(PCA_pred),
  function(x) unique(Use_cell_seurat@meta.data[Use_cell_seurat$orig.ident == x, "subtype"]),
  FUN.VALUE = character(1)
)

PCA_input <- data.frame(
  PC1 = PCA_pred[, 1],
  PC2 = PCA_pred[, 2],
  PC3 = PCA_pred[, 3],
  subtype = sample_subtype,
  row.names = rownames(PCA_pred),
  check.names = FALSE
)

PCA_input$color <- dplyr::case_when(
  PCA_input$subtype == "LME_1" ~ subtype_col1[1],
  PCA_input$subtype == "LME_2" ~ subtype_col1[2],
  PCA_input$subtype == "LME_3" ~ subtype_col1[3],
  TRUE ~ "grey80"
)

PCA_importance <- summary(PCA_prcomp)$importance["Proportion of Variance", ]

# ---- 3D scatter ----
pdf(file.path(fig_dir, "fig2e.pdf"), width = 8, height = 8)
scatter3D(
  x = PCA_input$PC1, y = PCA_input$PC3, z = PCA_input$PC2,
  xlab = paste0("PC1 (", round(PCA_importance[1], 3) * 100, "%)"),
  ylab = paste0("PC3 (", round(PCA_importance[3], 3) * 100, "%)"),
  zlab = paste0("PC2 (", round(PCA_importance[2], 3) * 100, "%)"),
  pch = 21, cex = 1.2,
  phi = 45, theta = 45,
  labels = PCA_input$subtype,
  col = "white",
  bg  = PCA_input$color
)
legend(
  "right",
  legend = c("LME 1", "LME 2", "LME 3"),
  col    = c(subtype_col1[1], subtype_col1[2], subtype_col1[3]),
  pch = 19, bty = "n"
)
dev.off()

# ---- build lightweight Seurat object with raw counts ----
counts_matrix <- Use_cell_seurat@assays$RNA@counts
new_use_cell_seurat <- CreateSeuratObject(counts = counts_matrix)

# carry over metadata (including fine_cell_type)
new_use_cell_seurat@meta.data <- Use_cell_seurat@meta.data

# ensure cell type labels are plain characters for export compatibility
new_use_cell_seurat$fine_cell_type <- as.character(new_use_cell_seurat$fine_cell_type)

# ---- export reference object ----
# First save as h5Seurat, then convert to h5ad format
SaveH5Seurat(new_use_cell_seurat,
             filename = "snrna_reference.h5Seurat",
             overwrite = TRUE)
Convert("snrna_reference.h5Seurat",
        dest = "h5ad",
        overwrite = TRUE)

# ---- note on intent ----
# The exported .h5ad file is specifically prepared as the single-cell reference
# for cell2location in spatial RNA-seq analysis.
