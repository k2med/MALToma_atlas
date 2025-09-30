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

# ---- derive CM1 / CM2 subsets from the merged object ----
cm1_members <- c("cDC_XCR1","pDC_CLEC4C","FDC_CR2","Mac_SPARC","Fibro_CCL19",
                 "cDC_CLEC10A","Endo_CCL21","CD4T_CXCL13","CD8T_HAVCR2","CD4T_FOXP3")
cm2_members <- c("Mac_MARCO","Endo_ACKR1","CD8T_GNLY","Mac_STAB1","CD8T_CCL5",
                 "Neu_G0S2","Endo_CXCL12","Peri_PDGFRB")

cellchat_cm1 <- subsetCellChat(cellchat_merged, idents.use = cm1_members)
cellchat_cm2 <- subsetCellChat(cellchat_merged, idents.use = cm2_members)

# ---- plots: interaction weights ----
p_total <- compareInteractions(
  cellchat_merged, show.legend = FALSE, group = c(1,2,3), measure = "weight",
  width = 0.5, size.text = 5
) +
  scale_fill_manual(values = unname(subtype_col1)) +
  theme(
    text = element_text(size = 5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.116),
    axis.ticks = element_line(colour = "black", linewidth = 0.116),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 5),
    plot.title = element_text(size = 5, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5)
  ) +
  coord_cartesian(clip = "off")

ggsave(
  filename = file.path(fig_dir, "fig6c_left.pdf"),
  plot     = egg::set_panel_size(p_total, width = unit(1, "cm"), height = unit(1.25, "cm")),
  width = 8, height = 5, units = "cm"
)

p_cm1 <- compareInteractions(
  cellchat_cm1, show.legend = FALSE, group = c(1,2,3), measure = "weight",
  width = 0.5, size.text = 5
) +
  scale_fill_manual(values = unname(subtype_col1)) +
  theme(
    text = element_text(size = 5),
    strip.text = element_text(size = 5),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.116),
    axis.ticks = element_line(colour = "black", linewidth = 0.116),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 5),
    plot.title = element_text(size = 5, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5)
  ) +
  coord_cartesian(clip = "off")

ggsave(
  filename = file.path(fig_dir, "fig6c_middle.pdf"),
  plot     = egg::set_panel_size(p_cm1, width = unit(1, "cm"), height = unit(1.25, "cm")),
  width = 8, height = 5, units = "cm"
)

p_cm2 <- compareInteractions(
  cellchat_cm2, show.legend = FALSE, group = c(1,2,3), measure = "weight",
  width = 0.5, size.text = 5
) +
  scale_fill_manual(values = unname(subtype_col1)) +
  theme(
    text = element_text(size = 5),
    strip.text = element_text(size = 5),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.116),
    axis.ticks = element_line(colour = "black", linewidth = 0.116),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 5),
    plot.title = element_text(size = 5, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5)
  ) +
  coord_cartesian(clip = "off")

ggsave(
  filename = file.path(fig_dir, "fig6c_right.pdf"),
  plot     = egg::set_panel_size(p_cm2, width = unit(1, "cm"), height = unit(1.25, "cm")),
  width = 8, height = 5, units = "cm"
)

# ---- pathway ranking comparison ----
use_pathway <- c("VEGF","TGFb","PDGF","COLLAGEN","CCL","ICAM","CXCL","CD40","BAFF")

p_rank <- rankNet(
  cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE,
  comparison = c(1,2,3), signaling = use_pathway
) +
  scale_x_discrete(limits = rev(use_pathway)) +
  scale_fill_manual(values = subtype_col2) +
  theme(
    text = element_text(size = 5),
    strip.text = element_text(size = 5),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.116),
    axis.ticks = element_line(colour = "black", linewidth = 0.116),
    axis.text  = element_text(colour = "black", size = 5),
    plot.title = element_text(size = 5),
    legend.position = "right",
    legend.title    = element_text(size = 5),
    legend.text     = element_text(size = 5)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.5),
    labels = function(x) x * 100
  )

ggsave(
  filename = file.path(fig_dir, "fig6d.pdf"),
  plot     = egg::set_panel_size(p_rank, width = unit(2, "cm"), height = unit(2, "cm")),
  width = 8, height = 5, units = "cm"
)
