# ---- packages ----
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(egg)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- inputs: main objects ----
seurat_myeloid <- readRDS(file.path(data_dir, "snrna_myeloid_subset_seurat.rds"))

# ---- subsets ----
seurat_mac <- subset(
  seurat_myeloid,
  fine_cell_type %in% c("Mac_SPARC", "Mac_STAB1", "Mac_MARCO", "Mac_CHIT1")
)
seurat_mac$fine_cell_type <- droplevels(seurat_mac$fine_cell_type)

# ---- UMAP on Harmony embedding ----
seurat_mac <- RunUMAP(seurat_mac, reduction = "harmony", dims = 1:20, min.dist = 0.3)

mac_umap <- DimPlot(
  object  = seurat_mac,
  pt.size = 2,
  group.by = "fine_cell_type",
  cols = fine_cell_type_col,
  raster = FALSE
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

ggsave(
  filename = file.path(fig_dir, "fig4g_left.pdf"),
  plot = egg::set_panel_size(
    mac_umap,
    width  = grid::unit(16, "cm"),
    height = grid::unit(16, "cm")
  ),
  width = 24, height = 24, units = "cm"
)

# ---- feature plots (MRC1 / FOLR2) ----
for (gene in c("MRC1", "FOLR2")) {
  feature_plot <- FeaturePlot(
    object  = seurat_mac,
    pt.size = 2,
    features = gene,
    cols = colorRampPalette(c("#dddddd", rev(hcl.colors(n = 10, palette = "Earth")[1:6])))(20),
    raster = FALSE
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
  
  ggsave(
    filename = file.path(fig_dir, paste0("fig4g_right_", gene, ".pdf")),
    plot = egg::set_panel_size(
      feature_plot,
      width  = grid::unit(16, "cm"),
      height = grid::unit(16, "cm")
    ),
    width = 24, height = 24, units = "cm"
  )
}

# ---- differential expression: CD206-axes & FOLR2-axes contrasts ----
cd206_markers <- FindMarkers(
  seurat_mac,
  ident.1 = c("Mac_STAB1", "Mac_MARCO"),
  ident.2 = c("Mac_SPARC", "Mac_CHIT1"),
  assay = "RNA",
  only.pos = FALSE,
  test.use = "MAST",
  min.pct = 0,
  logfc.threshold = 0
)

folr2_markers <- FindMarkers(
  seurat_mac,
  ident.1 = c("Mac_SPARC", "Mac_STAB1"),
  ident.2 = c("Mac_MARCO", "Mac_CHIT1"),
  assay = "RNA",
  only.pos = FALSE,
  test.use = "MAST",
  min.pct = 0,
  logfc.threshold = 0
)

colnames(cd206_markers) <- paste0("CD206_", colnames(cd206_markers))
colnames(folr2_markers) <- paste0("FOLR2_", colnames(folr2_markers))

merge_markers <- cbind(
  cd206_markers,
  folr2_markers[rownames(cd206_markers), ]
)

plot_input <- merge_markers[, c("CD206_avg_log2FC", "FOLR2_avg_log2FC")]
plot_input$gene <- rownames(plot_input)

# optional thresholds
logfc_cutoff <- 1
pct_cutoff   <- 0.5

gene_list <- c(
  "CXCL13", "SPARC", "ITGAX",
  "STAB1", "FOLR2",
  "MRC1", "CD163",
  "MARCO", "TREM2",
  "CHIT1", "LYZ"
)

macro_volcano <- ggplot(plot_input, aes(x = CD206_avg_log2FC, y = FOLR2_avg_log2FC)) +
  geom_point(size = 0.1, colour = "grey") +
  geom_point(
    data = subset(plot_input, gene %in% gene_list),
    size = 0.2, colour = "#bc1724"
  ) +
  theme_bw() +
  labs(x = "log2FC (MRC1+ vs. MRC1-)", y = "log2FC (FOLR2+ vs. FOLR2-)") +
  scale_x_continuous(limits = c(-4, 4)) +
  scale_y_continuous(limits = c(-4, 4)) +
  ggrepel::geom_text_repel(
    data = subset(plot_input, gene %in% gene_list),
    aes(label = gene),
    size = 5 / .pt,
    segment.size = 0.116,
    min.segment.length = 0.1,
    point.padding = 0.3,
    box.padding = 0.5,
    max.overlaps = 20,
    fontface = "italic"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = 3, colour = "black", linewidth = 0.116) +
  geom_hline(yintercept = c(-1, 1), linetype = 3, colour = "black", linewidth = 0.116) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.116),
    axis.line = element_line(colour = "black", linewidth = 0.116),
    axis.ticks = element_line(colour = "black", linewidth = 0.116),
    axis.text.x = element_text(colour = "black", size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    axis.title.x = element_text(colour = "black", size = 5),
    axis.title.y = element_text(colour = "black", size = 5),
    legend.title = element_text(size = 5)
  )

ggsave(
  filename = file.path(fig_dir, "fig4h.pdf"),
  plot = macro_volcano,
  width = 4, height = 4, units = "cm"
)