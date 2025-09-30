# ---- packages ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(CellChat)
library(egg)
library(grid)
library(future)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))
source(file.path(src_dir, "build_cellchat_merged.R"))
source(file.path(src_dir, "plot_signalingRole_heatmap_list.R"))

# ---- inputs ----
seurat_all <- readRDS(file.path(data_dir, "snrna_all_cells_seurat.rds"))

seurat_cm1_obj <- subset(
  seurat_all,
  fine_cell_type %in% c(
    "cDC_XCR1","pDC_CLEC4C","FDC_CR2","Mac_SPARC","Fibro_CCL19",
    "cDC_CLEC10A","Endo_CCL21","CD4T_CXCL13","CD8T_HAVCR2","CD4T_FOXP3"
  )
)
seurat_cm1_obj$coarse_cell_type <- droplevels(seurat_cm1_obj$coarse_cell_type)
seurat_cm1_obj$fine_cell_type   <- droplevels(seurat_cm1_obj$fine_cell_type)
Idents(seurat_cm1_obj) <- seurat_cm1_obj$fine_cell_type

seurat_cm2_obj <- subset(
  seurat_all,
  fine_cell_type %in% c(
    "Mac_MARCO","Endo_ACKR1","CD8T_GNLY","Mac_STAB1","CD8T_CCL5",
    "Neu_G0S2","Endo_CXCL12","Peri_PDGFRB"
  )
)
seurat_cm2_obj$coarse_cell_type <- droplevels(seurat_cm2_obj$coarse_cell_type)
seurat_cm2_obj$fine_cell_type   <- droplevels(seurat_cm2_obj$fine_cell_type)
Idents(seurat_cm2_obj) <- seurat_cm2_obj$fine_cell_type

# ---- run CellChat end-to-end in memory ----
cc_cm1_res <- build_cellchat_merged(
  seurat_obj     = seurat_cm1_obj,
  subtype_levels = c("LME_1","LME_2","LME_3"),
  group_col      = "fine_cell_type",
  subtype_col    = "subtype",
  assay          = "RNA",
  species        = "human",
  min_cells      = 10,
  workers        = 1
)
cellchat_cm1_merged <- cc_cm1_res$merged
cellchat_cm1_list   <- cc_cm1_res$per_subtype

cc_cm2_res <- build_cellchat_merged(
  seurat_obj     = seurat_cm2_obj,
  subtype_levels = c("LME_1","LME_2","LME_3"),
  group_col      = "fine_cell_type",
  subtype_col    = "subtype",
  assay          = "RNA",
  species        = "human",
  min_cells      = 10,
  workers        = 1
)
cellchat_cm2_merged <- cc_cm2_res$merged
cellchat_cm2_list   <- cc_cm2_res$per_subtype

# ---- helper: differential signaling-role scatter and save ----
save_diff_scatter <- function(cellchat_obj, comparison_vec, outfile) {
  p <- netAnalysis_diff_signalingRole_scatter(
    object     = cellchat_obj,
    comparison = comparison_vec,
    label.size = 5 / .pt,
    dot.size   = 0.5,
    do.label   = FALSE
  ) +
    scale_colour_manual(values = fine_cell_type_col, drop = FALSE) +
    theme(
      text            = element_text(size = 5),
      panel.grid      = element_blank(),
      panel.background= element_blank(),
      panel.border    = element_rect(colour = "black", linewidth = 0.116),
      axis.ticks      = element_line(colour = "black", linewidth = 0.116),
      axis.text.x     = element_text(colour = "black", size = 5),
      axis.text.y     = element_text(colour = "black", size = 5),
      plot.title      = element_text(size = 5, hjust = 0.5),
      legend.position = "right",
      legend.title    = element_text(size = 5),
      legend.text     = element_text(size = 5)
    ) +
    ggrepel::geom_text_repel(
      mapping       = aes(label = labels, colour = labels),
      size          = 5 / .pt,
      show.legend   = FALSE,
      segment.size  = 0.2,
      segment.alpha = 0.5,
      box.padding   = 0
    )
  
  ggsave(
    filename = file.path(fig_dir, outfile),
    plot     = egg::set_panel_size(p, width = unit(2, "cm"), height = unit(2, "cm")),
    width    = 5, height = 5, units = "cm"
  )
}

# ---- figures ----
# CM1: LME_2 vs LME_1 (top-left), LME_2 vs LME_3 (top-right)
save_diff_scatter(cellchat_cm1_merged, comparison_vec = c(1, 2), outfile = "fig6e_top_left.pdf")
save_diff_scatter(cellchat_cm1_merged, comparison_vec = c(3, 2), outfile = "fig6e_top_right.pdf")

# CM1: overall signaling-role heatmaps for CXCL / CD40 / ICAM / CCL
pdf(file.path(fig_dir, "fig6f.pdf"), width = 10, height = 4)
plot_signalingRole_heatmap_list(
  object.list   = cellchat_cm1_list,
  signaling     = c("CXCL", "CD40", "ICAM", "CCL"),
  pattern       = "all",
  width         = 2,
  height        = 0.8,
  color.heatmap = "YlOrRd",
  color.bar     = "#0571b0",
  color.use     = fine_cell_type_col[levels(cellchat_cm1_list$LME_1@idents)],
  font.size     = 5,
  font.size.title = 5
)
dev.off()

# CM2: LME_2 vs LME_1 (bottom-left), LME_2 vs LME_3 (bottom-right)
save_diff_scatter(cellchat_cm2_merged, comparison_vec = c(1, 2), outfile = "fig6e_bottom_left.pdf")
save_diff_scatter(cellchat_cm2_merged, comparison_vec = c(3, 2), outfile = "fig6e_bottom_right.pdf")

# CM2: overall signaling-role heatmaps for VEGF / TGF-β / PDGF / Collagen
pdf(file.path(fig_dir, "fig6h.pdf"), width = 10, height = 4)
plot_signalingRole_heatmap_list(
  object.list   = cellchat_cm2_list,
  signaling     = c("VEGF", "TGFb", "PDGF", "COLLAGEN"),
  pattern       = "all",
  width         = 1.6,
  height        = 0.8,
  color.heatmap = "YlOrRd",
  color.bar     = "#0571b0",
  color.use     = fine_cell_type_col[levels(cellchat_cm2_list$LME_1@idents)],
  font.size     = 5,
  font.size.title = 5
)
dev.off()