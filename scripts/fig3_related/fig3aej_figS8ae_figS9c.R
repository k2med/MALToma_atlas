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
seurat_stromal  <- readRDS(file.path(data_dir, "snrna_stromal_subset_seurat.rds"))

# ---- subsets ----
seurat_cd8t <- subset(seurat_t, fine_cell_type %in% c("CD8T_Teff", "CD8T_Temra", "CD8T_Tex")) %>%
  { .$fine_cell_type <- droplevels(.$fine_cell_type); . }

seurat_cd4t <- subset(seurat_t, fine_cell_type %in% c("CD4T_Treg", "CD4T_Tfh", "CD4T_Tn")) %>%
  { .$fine_cell_type <- droplevels(.$fine_cell_type); . }

seurat_mac <- subset(seurat_myeloid, fine_cell_type %in% c("Mac_SPARC", "Mac_STAB1", "Mac_MARCO", "Mac_CHIT1")) %>%
  { .$fine_cell_type <- droplevels(.$fine_cell_type); . }

seurat_dc <- subset(seurat_myeloid, fine_cell_type %in% c("cDC1_XCR1", "cDC2_CLEC10A", "pDC_CLEC4C")) %>%
  { .$fine_cell_type <- droplevels(.$fine_cell_type); . }

seurat_fibro <- subset(seurat_stromal, fine_cell_type %in% c("Fibro_PI16", "Fibro_CCL19", "Fibro_LRRC15")) %>%
  { .$fine_cell_type <- droplevels(.$fine_cell_type); . }

seurat_endo  <- subset(seurat_stromal, fine_cell_type %in% c("Endo_venous", "Endo_arterial", "Endo_lymphatic")) %>%
  { .$fine_cell_type <- droplevels(.$fine_cell_type); . }

# ---- outputs ----
plot_cell_type_proportion(
  seurat_obj       = seurat_cd8t,
  all_seurat_meta  = seurat_all@meta.data,
  cell_type_col    = "fine_cell_type",
  group_col        = "subtype",
  orig_col         = "orig.ident",
  cell_type_colors = fine_cell_type_col,
  filename         = "fig3a.pdf"
)

plot_cell_type_proportion(
  seurat_obj       = seurat_cd4t,
  all_seurat_meta  = seurat_all@meta.data,
  cell_type_col    = "fine_cell_type",
  group_col        = "subtype",
  orig_col         = "orig.ident",
  cell_type_colors = fine_cell_type_col,
  filename         = "figS8a.pdf"
)

plot_cell_type_proportion(
  seurat_obj       = seurat_mac,
  all_seurat_meta  = seurat_all@meta.data,
  cell_type_col    = "fine_cell_type",
  group_col        = "subtype",
  orig_col         = "orig.ident",
  cell_type_colors = fine_cell_type_col,
  filename         = "fig3e.pdf"
)

plot_cell_type_proportion(
  seurat_obj       = seurat_dc,
  all_seurat_meta  = seurat_all@meta.data,
  cell_type_col    = "fine_cell_type",
  group_col        = "subtype",
  orig_col         = "orig.ident",
  cell_type_colors = fine_cell_type_col,
  filename         = "figS8e.pdf"
)

plot_cell_type_proportion(
  seurat_obj       = seurat_fibro,
  all_seurat_meta  = seurat_all@meta.data,
  cell_type_col    = "fine_cell_type",
  group_col        = "subtype",
  orig_col         = "orig.ident",
  cell_type_colors = fine_cell_type_col,
  filename         = "fig3j.pdf"
)

plot_cell_type_proportion(
  seurat_obj       = seurat_endo,
  all_seurat_meta  = seurat_all@meta.data,
  cell_type_col    = "fine_cell_type",
  group_col        = "subtype",
  orig_col         = "orig.ident",
  cell_type_colors = fine_cell_type_col,
  filename         = "figS9c.pdf"
)