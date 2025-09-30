# ---- packages ----
library(Seurat)
library(AUCell)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(egg)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- inputs: Seurat objects ----
seurat_t <- readRDS(file.path(data_dir, "snrna_t_subset_seurat.rds"))
seurat_cd4t <- subset(
  seurat_t,
  fine_cell_type %in% c("CD4T_FOXP3", "CD4T_CXCL13", "CD4T_IL7R")
)

# ---- compute Treg proportion per sample ----
treg_prop_df <- data.frame(
  Treg_prop = prop.table(
    table(seurat_cd4t$orig.ident, seurat_cd4t$fine_cell_type == "CD4T_FOXP3"),
    margin = 1
  )[, 2]
)

# ---- inputs: B_Mal Seurat and AUCell scores from fig3cd.R ----
# cells_auc was generated in fig3cd.R and saved as "cells_AUC.rds"
cells_auc <- readRDS(file.path(fig_dir, "cells_AUC.rds"))

seurat_b <- readRDS(file.path(data_dir, "snrna_b_subset_seurat.rds"))
seurat_bmal <- subset(seurat_b, fine_cell_type == "B_Mal")

auc_mat <- AUCell::getAUC(cells_auc)  # genesets x cells
stopifnot(identical(rownames(seurat_bmal@meta.data), colnames(auc_mat)))

meta_df <- cbind(seurat_bmal@meta.data, t(auc_mat))

# ---- aggregate MP scores per sample (median) ----
mp_median_by_sample <- meta_df %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(dplyr::across(MP_1:MP_7, median, na.rm = FALSE)) %>%
  tibble::column_to_rownames(var = "orig.ident")

# align and merge with Treg proportion
treg_prop_mp_df <- cbind(treg_prop_df, mp_median_by_sample[rownames(treg_prop_df), , drop = FALSE])

# add subtype per sample (from meta_df)
treg_prop_mp_df$subtype <- sapply(
  rownames(treg_prop_mp_df),
  function(x) unique(meta_df[meta_df$orig.ident == x, "subtype"])
)

# ---- plot: correlation scatter (MP_7 vs Treg proportion) ----
cor_scatter <- ggplot(treg_prop_mp_df, aes(x = MP_7, y = Treg_prop, color = subtype)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.116, color = "black") +
  stat_cor(
    method  = "pearson",
    label.x = min(treg_prop_mp_df$MP_7, na.rm = TRUE),
    label.y = max(treg_prop_mp_df$Treg_prop, na.rm = TRUE),
    size    = 5 / .pt,
    color   = "black"
  ) +
  scale_color_manual(values = subtype_col1) +
  theme_bw(base_size = 5) +
  labs(x = "MP_7 score", y = "Proportion of CD4T_FOXP3") +
  theme(
    text = element_text(size = 5),
    panel.grid.major = element_line(colour = "grey90", size = 0.116),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", size = 0.116),
    axis.line = element_line(colour = "black", size = 0.116),
    axis.ticks = element_line(colour = "black", size = 0.116),
    axis.text.x = element_text(colour = "black", size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    legend.position = "right",
    legend.title = element_text(size = 5)
  )

ggsave(
  filename = file.path(fig_dir, "fig3g.pdf"),
  plot = egg::set_panel_size(
    cor_scatter,
    width  = grid::unit(2, "cm"),
    height = grid::unit(2, "cm")
  ),
  width = 10, height = 10, units = "cm"
)
