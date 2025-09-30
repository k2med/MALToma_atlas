# ---- packages ----
library(beyondcell)
library(Seurat)
library(tidyverse)
library(ggridges)
library(egg)
library(grid)

# ---- paths & sources ----
snrna_data_dir <- "../../data/snrna_rnaseq"
src_dir        <- "../source"
fig_dir        <- "."

source(file.path(src_dir, "custom_colors.R"))
source(file.path(src_dir, "beyondcell_bcRanks.R"))
source(file.path(src_dir, "plot_beyondcell_ridgeplot.R"))

# ---- inputs ----
seurat_b    <- readRDS(file.path(snrna_data_dir, "snrna_b_subset_seurat.rds"))
seurat_bmal <- subset(seurat_b, fine_cell_type == "B_Mal")

# ---- beyondcell scoring & ranking ----
gene_sets    <- GetCollection(SSc)
bc_scores    <- bcScore(seurat_bmal, gene_sets, expr.thres = 0.1)
bc_scores    <- bcRanks(bc_scores, idents = "subtype")

# ---- ridge plots helper ----
save_ridge <- function(sig_id, outfile) {
  ridge_plot <- bcRidge(bc_scores, signatures = sig_id, idents = "subtype") +
    scale_fill_manual(values = subtype_col1)
  ggsave(
    filename = file.path(fig_dir, outfile),
    plot = egg::set_panel_size(
      ridge_plot,
      width  = grid::unit(1.5, "cm"),
      height = grid::unit(1.5, "cm")
    ),
    width = 12, height = 6, units = "cm"
  )
}

# ---- outputs ----
save_ridge("sig-21115", "fig3i_sorafenib.pdf")
save_ridge("sig-21303", "fig3i_ibrutinib.pdf")
save_ridge("sig-21353", "fig3i_vincristine.pdf")
save_ridge("sig-21139", "figS3d_axitinib.pdf")
save_ridge("sig-21175", "figS3d_sunitinib.pdf")
save_ridge("sig-20958", "figS3d_vinorelbine.pdf")
save_ridge("sig-21025", "figS3d_tak901.pdf")