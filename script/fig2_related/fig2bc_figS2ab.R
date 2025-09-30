# ---- packages ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(egg)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- all cells ----
All_cell_seurat <- readRDS(file.path(data_dir, "snrna_all_cells_seurat.rds"))

p_all_umap <- DimPlot(
  object   = All_cell_seurat,
  pt.size  = 0.1,
  group.by = "coarse_cell_type",
  cols     = coarse_cell_type_col,
  raster   = FALSE
) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    strip.text   = element_blank(),
    aspect.ratio = 1
  )

ggsave(file.path(fig_dir, "fig2b.pdf"), plot = p_all_umap,
       width = 18, height = 16, units = "cm")

# ---- coarse marker dotplot ----
Idents(All_cell_seurat) <- All_cell_seurat$coarse_cell_type
coarse_markers <- c("CD79A","MS4A1","TRAC","TRBC1",
                    "CD68","LYZ","COL1A1","VWF",
                    "ALB","APOA1","AGER","SFTPB")

dp <- DotPlot(All_cell_seurat, features = coarse_markers)
dotplot_data <- dp$data[, c("id","features.plot","pct.exp","avg.exp.scaled")]
dotplot_data$id <- factor(dotplot_data$id, levels = rev(levels(dotplot_data$id)))

coarse_dotplot <- ggplot(dotplot_data, aes(x = id, y = features.plot)) +
  geom_point(aes(color = avg.exp.scaled, size = pct.exp),
             shape = 19, stroke = 0.01) +
  xlab("") + ylab("") +
  scale_color_gradientn(colors = colorRampPalette(c("#ffffff","#fee0d2","#ae0000"))(20)) +
  scale_size(range = c(0, 3.2), limits = c(0, 100),
             breaks = c(0,20,40,60,80,100)) +
  theme_bw(base_size = 5) +
  theme(
    text         = element_text(size = 5),
    panel.grid   = element_blank(),
    panel.border = element_rect(colour = "black", size = 0.116),
    axis.line    = element_line(colour = "black", size = 0.116),
    axis.ticks   = element_line(colour = "black", size = 0.116),
    axis.text.x  = element_text(colour = "black", angle = 60, vjust = 1, hjust = 1,
                                size = 5, face = "italic"),
    axis.text.y  = element_text(colour = "black", size = 5),
    legend.position = "right",
    legend.title    = element_text(size = 5)
  ) +
  guides(size = guide_legend(ncol = 1, byrow = TRUE,
                             override.aes = list(stroke = 0.4)),
         fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  coord_flip()

ggsave(file.path(fig_dir, "figS2a.pdf"),
       egg::set_panel_size(coarse_dotplot,
                           width  = unit(length(coarse_markers)/3, "cm"),
                           height = unit(length(levels(dotplot_data$id))/3, "cm")),
       width = 10, height = 10, units = "cm")

# ---- B cells ----
B_cell_seurat <- readRDS(file.path(data_dir, "snrna_b_subset_seurat.rds"))
B_cell_seurat$fine_cell_type <- factor(B_cell_seurat$fine_cell_type,
                                       levels = c("B_Mal","B_Nor"))

B_UMAP <- DimPlot(B_cell_seurat, pt.size = 0.1,
                  group.by = "fine_cell_type",
                  cols = fine_cell_type_col, raster = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    strip.text   = element_blank(),
    aspect.ratio = 1
  )

ggsave(file.path(fig_dir, "fig2c_b_cell.pdf"),
       egg::set_panel_size(B_UMAP, width = unit(16,"cm"), height = unit(16,"cm")),
       width = 24, height = 24, units = "cm")

# ---- T cells ----
T_cell_seurat <- readRDS(file.path(data_dir, "snrna_t_subset_seurat.rds"))
T_cell_seurat$fine_cell_type <- factor(T_cell_seurat$fine_cell_type,
                                       levels = c("CD8T_CCL5","CD8T_GNLY","CD8T_HAVCR2",
                                                  "CD4T_FOXP3","CD4T_CXCL13","CD4T_IL7R"))

T_UMAP <- DimPlot(T_cell_seurat, pt.size = 0.5,
                  group.by = "fine_cell_type",
                  cols = fine_cell_type_col, raster = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    strip.text   = element_blank(),
    aspect.ratio = 1
  )

ggsave(file.path(fig_dir, "fig2c_t_cell.pdf"),
       egg::set_panel_size(T_UMAP, width = unit(16,"cm"), height = unit(16,"cm")),
       width = 24, height = 24, units = "cm")

# ---- Myeloid ----
Myeloid_cell_seurat <- readRDS(file.path(data_dir, "snrna_myeloid_subset_seurat.rds"))
Myeloid_cell_seurat$fine_cell_type <- factor(Myeloid_cell_seurat$fine_cell_type,
                                             levels = c("Mac_SPARC","Mac_STAB1","Mac_MARCO","Mac_CHIT1",
                                                        "cDC_XCR1","cDC_CLEC10A","pDC_CLEC4C",
                                                        "Neu_G0S2","Mast_CPA3"))

Myeloid_UMAP <- DimPlot(Myeloid_cell_seurat, pt.size = 0.5,
                        group.by = "fine_cell_type",
                        cols = fine_cell_type_col, raster = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    strip.text   = element_blank(),
    aspect.ratio = 1
  )

ggsave(file.path(fig_dir, "fig2c_myeloid_cell.pdf"),
       egg::set_panel_size(Myeloid_UMAP, width = unit(16,"cm"), height = unit(16,"cm")),
       width = 24, height = 24, units = "cm")

# ---- Stromal ----
Stromal_cell_seurat <- readRDS(file.path(data_dir, "snrna_stromal_subset_seurat.rds"))
Stromal_cell_seurat$fine_cell_type <- factor(Stromal_cell_seurat$fine_cell_type,
                                             levels = c("Endo_ACKR1","Endo_CXCL12","Endo_CCL21","Endo_HPGD","Peri_PDGFRB",
                                                        "Fibro_PI16","Fibro_CCL19","Fibro_LRRC15","Fibro_TCF21","FDC_CR2"))

Stromal_UMAP <- DimPlot(Stromal_cell_seurat, pt.size = 0.5,
                        group.by = "fine_cell_type",
                        cols = fine_cell_type_col, raster = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    strip.text   = element_blank(),
    aspect.ratio = 1
  )

ggsave(file.path(fig_dir, "fig2c_stromal_cell.pdf"),
       egg::set_panel_size(Stromal_UMAP, width = unit(16,"cm"), height = unit(16,"cm")),
       width = 24, height = 24, units = "cm")

# ---- Alveolar ----
Alveolar_cell_seurat <- readRDS(file.path(data_dir, "snrna_alveolar_subset_seurat.rds"))
Alveolar_cell_seurat$fine_cell_type <- factor(Alveolar_cell_seurat$fine_cell_type,
                                              levels = c("AT1_AGER","AT2_SFTPC"))

Alveolar_UMAP <- DimPlot(Alveolar_cell_seurat, pt.size = 1,
                         group.by = "fine_cell_type",
                         cols = fine_cell_type_col, raster = FALSE) +
  labs(x = "UMAP 1", y = "UMAP 2", title = NULL) +
  theme(
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    strip.text   = element_blank(),
    aspect.ratio = 1
  )

ggsave(file.path(fig_dir, "figS2b.pdf"),
       egg::set_panel_size(Alveolar_UMAP, width = unit(16,"cm"), height = unit(16,"cm")),
       width = 24, height = 24, units = "cm")