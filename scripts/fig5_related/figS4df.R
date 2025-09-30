# ---- packages ----
library(Seurat)
library(tidyverse)
library(CytoTRACE2)
library(ggplot2)
library(ggpubr)
library(egg)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- inputs ----
seurat_stromal <- readRDS(file.path(data_dir, "snrna_stromal_subset_seurat.rds"))

# ---- subsets ----
seurat_fibro <- subset(seurat_stromal, fine_cell_type %in% c("Fibro_PI16", "Fibro_CCL19", "Fibro_LRRC15"))
seurat_vec   <- subset(seurat_stromal, fine_cell_type %in% c("Endo_ACKR1", "Endo_CXCL12"))
seurat_fibro$fine_cell_type <- droplevels(seurat_fibro$fine_cell_type)
seurat_vec$fine_cell_type   <- droplevels(seurat_vec$fine_cell_type)

# ---- helper: run CytoTRACE2 + boxplot ----
run_ct2_boxplot <- function(seurat_obj, outfile, seed = 421) {
  ct2 <- cytotrace2(
    seurat_obj,
    is_seurat = TRUE,
    slot_type = "counts",
    species   = "human",
    seed      = seed
  )
  
  df <- ct2@meta.data[, c("fine_cell_type", "CytoTRACE2_Score")]
  p  <- ggplot(df, aes(x = fine_cell_type, y = CytoTRACE2_Score, fill = fine_cell_type)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.116, width = 0.6) +
    scale_fill_manual(values = fine_cell_type_col) +
    theme_minimal(base_size = 5) +
    theme(
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      axis.line        = element_line(colour = "black", linewidth = 0.116),
      axis.ticks       = element_line(colour = "black", linewidth = 0.116),
      axis.text.x      = element_text(colour = "black", size = 5, angle = 60, vjust = 1, hjust = 1),
      axis.ticks.x     = element_blank(),
      axis.text.y      = element_text(colour = "black", size = 5),
      legend.position  = "none",
      plot.title       = element_text(size = 5, hjust = 0.5)
    ) +
    ggpubr::stat_pwc(
      aes(group = fine_cell_type),
      label.size       = 5 * 25.4 / 72,
      tip.length       = 0,
      size             = 0.116,
      hide.ns          = TRUE,
      vjust            = 0.5,
      method           = "wilcox.test",
      p.adjust.method  = "BH",
      label            = "p.adj.signif"
    ) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "CytoTRACE2 score")
  
  ggsave(
    filename = file.path(fig_dir, outfile),
    plot     = egg::set_panel_size(p, width = grid::unit(1.5, "cm"), height = grid::unit(2, "cm")),
    width    = 8, height = 5, units = "cm"
  )
  invisible(p)
}

# ---- runs ----
run_ct2_boxplot(seurat_fibro, outfile = "figS4d.pdf")
run_ct2_boxplot(seurat_vec,   outfile = "figS4f.pdf")
