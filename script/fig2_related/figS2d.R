# ---- packages ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(egg)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))  # expects site_col

# ---- all cells ----
All_cell_seurat <- readRDS(file.path(data_dir, "snrna_all_cells_seurat.rds"))
meta_df <- All_cell_seurat@meta.data

# rows = fine_cell_type, cols = site; column-wise proportions then transpose
celltype_site_prop <- t(round(prop.table(table(meta_df$site, meta_df$fine_cell_type), margin = 2), 3))

# ---- Shannon entropy per cell type ----
entropy_df <- data.frame(
  shannon_entropy = apply(celltype_site_prop, 1, function(p) {
    p <- p[p > 0]
    -sum(p * log(p))
  })
)
entropy_df <- entropy_df[order(entropy_df$shannon_entropy, decreasing = TRUE), , drop = FALSE]
entropy_df$cell <- rownames(entropy_df)
entropy_df$cell <- factor(entropy_df$cell, levels = entropy_df$cell)
entropy_df$max_prop <- apply(celltype_site_prop[as.character(entropy_df$cell), , drop = FALSE], 1, max)

# ---- dot plot ----
entropy_dotplot <- ggplot(entropy_df, aes(x = cell, y = shannon_entropy)) +
  geom_point(aes(size = max_prop), color = "#517595") +
  scale_size(range = c(0, 2), limits = c(0.3, 1)) +
  labs(x = "Cell type", y = "Shannon entropy") +
  theme(
    text = element_text(size = 5),
    panel.grid.major.x = element_line(colour = "grey90", size = 0.116),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black", size = 0.116),
    axis.ticks.x       = element_blank(),
    axis.ticks.y       = element_line(colour = "black", size = 0.116),
    axis.title         = element_text(colour = "black", size = 5),
    axis.text.x        = element_blank(),
    axis.text.y        = element_text(colour = "black", size = 5),
    legend.position    = "bottom",
    legend.title       = element_text(size = 5),
    legend.text        = element_text(size = 5)
  ) +
  geom_text_repel(
    data = subset(entropy_df, shannon_entropy < 0.5),
    aes(label = cell),
    size = 5 * 25.4 / 72,
    segment.size = 0.116,
    min.segment.length = 0.1,
    color = "black",
    point.padding = 0.3
  )

ggsave(
  file.path(fig_dir, "figS2d_dotplot.pdf"),
  egg::set_panel_size(entropy_dotplot, width = unit(4, "cm"), height = unit(2.5, "cm")),
  width = 8, height = 6, units = "cm"
)

# ---- cell proportions ----
celltype_site_prop <- as.data.frame(prop.table(table(meta_df$site, meta_df$fine_cell_type), margin = 2))
colnames(celltype_site_prop) <- c("site", "cell", "freq")

# ---- tidy for pie ----
cells_keep <- c("Fibro_TCF21", "AT1_AGER", "AT2_SFTPC", "Hepa_ALB", "Endo_HPGD")

pie_df <- celltype_site_prop %>%
  filter(cell %in% cells_keep)

pie_df$site <- factor(
  pie_df$site,
  levels = c("Liver", "Lung", "Ocular_adnexa", "Paranasal_sinus", "Salivary_gland", "Lymph_node")
)
pie_df$cell <- factor(
  pie_df$cell,
  levels = c("Hepa_ALB", "AT2_SFTPC", "AT1_AGER", "Endo_HPGD", "Fibro_TCF21")
)

plot_pie <- ggplot(pie_df, aes(x = "", y = freq, fill = site)) +
  facet_grid(. ~ cell) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = site_col) +
  coord_polar(theta = "y") +
  theme_void(base_size = 5) +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 5),
    legend.text     = element_text(size = 5)
  )

ggsave(
  file.path(fig_dir, "figS2d_pie.pdf"),
  egg::set_panel_size(plot_pie, width = unit(1, "cm"), height = unit(1, "cm")),
  width = 10, height = 6, units = "cm"
)