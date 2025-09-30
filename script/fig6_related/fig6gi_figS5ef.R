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

# ---- inputs ----
seurat_all <- readRDS(file.path(data_dir, "snrna_all_cells_seurat.rds"))
seurat_use <- subset(
  seurat_all,
  fine_cell_type %in% c("Fibro_TCF21", "AT1_AGER", "AT2_SFTPC", "Hepa_ALB", "Endo_HPGD"),
  invert = TRUE
)
seurat_use$coarse_cell_type <- droplevels(seurat_use$coarse_cell_type)
seurat_use$fine_cell_type   <- droplevels(seurat_use$fine_cell_type)
Idents(seurat_use) <- seurat_use$fine_cell_type

# ---- run CellChat end-to-end in memory ----
cc_res <- build_cellchat_merged(
  seurat_obj     = seurat_use,
  subtype_levels = c("LME_1","LME_2","LME_3"),
  group_col      = "fine_cell_type",
  subtype_col    = "subtype",
  assay          = "RNA",
  species        = "human",
  min_cells      = 10,
  workers        = 1
)
cellchat_merged <- cc_res$merged

# ---- helper: bubble plot with consistent theming & saving ----
save_bubble_plot <- function(cellchat_obj,
                             sources, targets, signaling,
                             outfile,
                             comparison = c(1, 2, 3),
                             sort_source = TRUE,
                             sort_target = TRUE,
                             dot_min = 0, dot_max = 2,
                             angle_x = 45,
                             reverse_y = TRUE,
                             fig_w_cm = 16, fig_h_cm = 16,
                             panel_w_cm = NULL, panel_h_cm = NULL,
                             thresh = 0.05) {
  p <- netVisual_bubble(
    object          = cellchat_obj,
    sources.use     = sources,
    targets.use     = targets,
    signaling       = signaling,
    comparison      = comparison,
    sort.by.source  = sort_source,
    sort.by.target  = sort_target,
    dot.size.min    = dot_min,
    dot.size.max    = dot_max,
    angle.x         = angle_x,
    thresh          = thresh
  )
  
  p2 <- p +
    {
      if (reverse_y) {
        scale_y_discrete(limits = rev(levels(p$data$interaction_name_2)))
      } else {
        scale_y_discrete(limits = levels(p$data$interaction_name_2))
      }
    } +
    theme(
      text            = element_text(size = 5),
      strip.text      = element_text(size = 5),
      strip.background= element_blank(),
      panel.grid.major= element_blank(),
      panel.grid.minor= element_blank(),
      panel.background= element_blank(),
      axis.line       = element_line(colour = "black", linewidth = 0.116),
      axis.ticks      = element_line(colour = "black", linewidth = 0.116),
      axis.text       = element_text(colour = "black", size = 5),
      plot.title      = element_text(size = 5, hjust = 0.5),
      legend.position = "right",
      legend.title    = element_text(size = 5),
      legend.text     = element_text(size = 5)
    )
  
  # optional panel size for compact pdfs
  if (is.null(panel_w_cm)) panel_w_cm <- fig_w_cm
  if (is.null(panel_h_cm)) panel_h_cm <- fig_h_cm
  
  ggsave(
    filename = file.path(fig_dir, outfile),
    plot     = egg::set_panel_size(
      p2,
      width  = grid::unit(panel_w_cm, "cm"),
      height = grid::unit(panel_h_cm, "cm")
    ),
    width  = fig_w_cm,
    height = fig_h_cm,
    units  = "cm"
  )
}

# ---- Fig 6g: CXCL / ICAM ----
save_bubble_plot(
  cellchat_obj = cellchat_merged,
  sources  = c("Mac_SPARC", "FDC_CR2", "Fibro_CCL19"),
  targets  = c("CD4T_CXCL13", "CD4T_FOXP3", "CD8T_HAVCR2"),
  signaling= c("CXCL", "ICAM"),
  outfile  = "fig6g.pdf",
  dot_min  = 0, dot_max = 2,
  reverse_y   = TRUE,
  fig_w_cm    = 16, fig_h_cm = 16,
  panel_w_cm  = 24 * 0.25,
  panel_h_cm  = 7 * 0.25,
  thresh      = 0.05
)

# ---- Fig S5e: CD40 ----
save_bubble_plot(
  cellchat_obj = cellchat_merged,
  sources  = c("Mac_SPARC", "FDC_CR2", "Fibro_CCL19", "cDC_CLEC10A"),
  targets  = c("CD4T_CXCL13", "CD4T_FOXP3", "CD8T_HAVCR2"),
  signaling= "CD40",
  outfile  = "figS5e.pdf",
  dot_min  = 0, dot_max = 2,
  reverse_y   = FALSE,
  fig_w_cm    = 20, fig_h_cm = 16,
  panel_w_cm  = 36 * 0.25,
  panel_h_cm  = 3 * 0.25,
  thresh      = 0.05
)

# ---- Fig 6i: VEGF ----
save_bubble_plot(
  cellchat_obj = cellchat_merged,
  sources  = c("Neu_G0S2", "Mac_MARCO"),
  targets  = c("Endo_ACKR1", "Endo_CXCL12"),
  signaling= "VEGF",
  outfile  = "fig6i.pdf",
  dot_min  = 1, dot_max = 2,
  reverse_y   = TRUE,
  fig_w_cm    = 16, fig_h_cm = 16,
  panel_w_cm  = 12 * 0.25,
  panel_h_cm  = 9 * 0.25,
  thresh      = 0.05
)

# ---- Fig S5f: TGFb ----
save_bubble_plot(
  cellchat_obj = cellchat_merged,
  sources  = c("Neu_G0S2", "Mac_MARCO"),
  targets  = c("Endo_ACKR1", "Endo_CXCL12"),
  signaling= "TGFb",
  outfile  = "figS5f.pdf",
  dot_min  = 1, dot_max = 2,
  reverse_y   = TRUE,
  fig_w_cm    = 16, fig_h_cm = 16,
  panel_w_cm  = 12 * 0.25,
  panel_h_cm  = 9 * 0.25,
  thresh      = 0.05
)