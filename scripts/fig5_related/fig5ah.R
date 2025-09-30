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
seurat_stromal  <- readRDS(file.path(data_dir, "snrna_stromal_subset_seurat.rds"))

# ---- subsets ----
seurat_fibro <- subset(seurat_stromal, fine_cell_type %in% c("Fibro_PI16", "Fibro_CCL19", "Fibro_LRRC15")) %>%
  { .$fine_cell_type <- droplevels(.$fine_cell_type); . }

seurat_endo  <- subset(seurat_stromal, fine_cell_type %in% c("Endo_ACKR1", "Endo_CXCL12", "Endo_CCL21")) %>%
  { .$fine_cell_type <- droplevels(.$fine_cell_type); . }

# ---- outputs ----
plot_cell_type_proportion(
  seurat_obj       = seurat_fibro,
  all_seurat_meta  = seurat_all@meta.data,
  cell_type_col    = "fine_cell_type",
  group_col        = "subtype",
  orig_col         = "orig.ident",
  cell_type_colors = fine_cell_type_col,
  filename         = "fig5a.pdf"
)

plot_cell_type_proportion(
  seurat_obj       = seurat_endo,
  all_seurat_meta  = seurat_all@meta.data,
  cell_type_col    = "fine_cell_type",
  group_col        = "subtype",
  orig_col         = "orig.ident",
  cell_type_colors = fine_cell_type_col,
  filename         = "fig5h.pdf"
)
