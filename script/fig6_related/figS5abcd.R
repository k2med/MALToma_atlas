# ---- packages ----
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(CellChat)
library(ComplexHeatmap)
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
object_list <- cc_res$per_subtype

# ---- line plot: BAFF / CD40 indegree on B_Mal across LME ----
lme_order <- c("LME_1", "LME_2", "LME_3")
pathways  <- c("CD40", "BAFF")

lme_values <- sapply(lme_order, function(lme) {
  c(
    CD40 = ifelse(is.null(object_list[[lme]]@netP$centr$CD40$indeg["B_Mal"]),
                  0, object_list[[lme]]@netP$centr$CD40$indeg["B_Mal"]),
    BAFF = ifelse(is.null(object_list[[lme]]@netP$centr$BAFF$indeg["B_Mal"]),
                  0, object_list[[lme]]@netP$centr$BAFF$indeg["B_Mal"])
  )
})

line_df <- data.frame(
  LME     = rep(lme_order, each = length(pathways)),
  Pathway = rep(pathways, times = length(lme_order)),
  Value   = as.vector(lme_values)
)

line_df$rel_strength <- -1 / log(line_df$Value)

line_plot <- ggplot(line_df, aes(x = LME, y = rel_strength, color = Pathway, group = Pathway)) +
  geom_line(linewidth = 0.116) +
  geom_point(size = 1, stroke = 0.116) +
  theme_minimal(base_size = 5) +
  scale_color_manual(values = c("#6f9bb9", "#f6cf78")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.116),
    axis.ticks       = element_line(colour = "black", linewidth = 0.116),
    axis.text.x      = element_text(colour = "black", size = 5, angle = 60, vjust = 1, hjust = 1),
    axis.text.y      = element_text(colour = "black", size = 5),
    legend.position  = "right",
    legend.title     = element_blank()
  ) +
  labs(x = NULL, y = "Relative strength", title = NULL) +
  coord_cartesian(clip = "off")

ggsave(
  filename = file.path(fig_dir, "figS5a.pdf"),
  plot     = egg::set_panel_size(line_plot, width = unit(2, "cm"), height = unit(1, "cm")),
  width = 10, height = 10, units = "cm"
)

# ---- outgoing signaling-role heatmaps for BAFF / CD40 ----
pdf(file.path(fig_dir, "figS5b.pdf"), width = 10, height = 4)
plot_signalingRole_heatmap_list(
  object.list   = object_list,
  signaling     = c("BAFF", "CD40"),
  pattern       = "outgoing",
  width         = 5,
  height        = 0.4,
  color.heatmap = "YlOrRd",
  color.bar     = "#0571b0",
  color.use     = fine_cell_type_col[levels(object_list$LME_1@idents)],
  font.size     = 5,
  font.size.title = 5
)
dev.off()

# ---- circle networks targeting B_Mal: BAFF then CD40 ----
plot_signal_circle <- function(obj_list, pathway, outfile, width = 6, height = 4) {
  weight_max <- getMaxWeight(obj_list, slot.name = "netP", attribute = pathway)
  pdf(file.path(fig_dir, outfile), width = width, height = height)
  par(mfrow = c(1, length(obj_list)), xpd = TRUE)
  for (i in seq_along(obj_list)) {
    netVisual_aggregate(
      object          = obj_list[[i]],
      signaling       = pathway,
      targets.use     = "B_Mal",
      layout          = "circle",
      remove.isolate  = TRUE,
      thresh          = 0.05,
      vertex.label.cex = 2/3,
      color.use       = fine_cell_type_col,
      edge.weight.max = weight_max[1],
      signaling.name  = paste0(pathway, " ", names(obj_list)[i])
    )
  }
  dev.off()
}

plot_signal_circle(object_list, pathway = "BAFF", outfile = "figS5c.pdf")
plot_signal_circle(object_list, pathway = "CD40", outfile = "figS5d.pdf")