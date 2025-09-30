# ---- packages ----
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggalluvial)
library(patchwork)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))
source(file.path(src_dir, "plot_cell_type_prop.R"))

# ---- inputs: main objects ----
seurat_all      <- readRDS(file.path(data_dir, "snrna_all_cells_seurat.rds"))
seurat_t        <- readRDS(file.path(data_dir, "snrna_t_subset_seurat.rds"))
seurat_myeloid  <- readRDS(file.path(data_dir, "snrna_myeloid_subset_seurat.rds"))

# ---- subsets ----
seurat_cd8t <- subset(seurat_t, fine_cell_type %in% c("CD8T_CCL5", "CD8T_GNLY", "CD8T_HAVCR2")) %>%
  { .$fine_cell_type <- droplevels(.$fine_cell_type); . }

seurat_cd4t <- subset(seurat_t, fine_cell_type %in% c("CD4T_FOXP3", "CD4T_CXCL13", "CD4T_IL7R")) %>%
  { .$fine_cell_type <- droplevels(.$fine_cell_type); . }

seurat_mac <- subset(seurat_myeloid, fine_cell_type %in% c("Mac_SPARC", "Mac_STAB1", "Mac_MARCO", "Mac_CHIT1")) %>%
  { .$fine_cell_type <- droplevels(.$fine_cell_type); . }

seurat_dc <- subset(seurat_myeloid, fine_cell_type %in% c("cDC_XCR1", "cDC_CLEC10A", "pDC_CLEC4C")) %>%
  { .$fine_cell_type <- droplevels(.$fine_cell_type); . }

# ---- outputs ----
plot_cell_type_proportion(
  seurat_obj       = seurat_cd8t,
  all_seurat_meta  = seurat_all@meta.data,
  cell_type_col    = "fine_cell_type",
  group_col        = "subtype",
  orig_col         = "orig.ident",
  cell_type_colors = fine_cell_type_col,
  filename         = "fig4a.pdf"
)

plot_cell_type_proportion(
  seurat_obj       = seurat_cd4t,
  all_seurat_meta  = seurat_all@meta.data,
  cell_type_col    = "fine_cell_type",
  group_col        = "subtype",
  orig_col         = "orig.ident",
  cell_type_colors = fine_cell_type_col,
  filename         = "fig4b.pdf"
)

plot_cell_type_proportion(
  seurat_obj       = seurat_mac,
  all_seurat_meta  = seurat_all@meta.data,
  cell_type_col    = "fine_cell_type",
  group_col        = "subtype",
  orig_col         = "orig.ident",
  cell_type_colors = fine_cell_type_col,
  filename         = "fig4f.pdf"
)

plot_cell_type_proportion(
  seurat_obj       = seurat_dc,
  all_seurat_meta  = seurat_all@meta.data,
  cell_type_col    = "fine_cell_type",
  group_col        = "subtype",
  orig_col         = "orig.ident",
  cell_type_colors = fine_cell_type_col,
  filename         = "fig4j.pdf"
)