# ---- packages ----
library(Seurat)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(egg)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- inputs ----
seurat_b <- readRDS(file.path(data_dir, "snrna_b_subset_seurat.rds"))
seurat_bmal <- subset(seurat_b, fine_cell_type == "B_Mal")

# ---- cell cycle scoring ----
s_genes   <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

seurat_bmal <- CellCycleScoring(
  seurat_bmal,
  s.features   = s_genes,
  g2m.features = g2m_genes,
  set.ident    = TRUE
)

# ---- plot: G2M score by subtype ----
g2m_score_boxplot <- ggplot(seurat_bmal@meta.data, aes(x = subtype, y = G2M.Score, fill = subtype)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.116, width = 0.6) +
  scale_fill_manual(values = subtype_col1) +
  coord_cartesian(ylim = c(-0.25, 0.2)) +
  theme_minimal(base_size = 5) +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.116),
    axis.ticks       = element_line(colour = "black", linewidth = 0.116),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(colour = "black", size = 5),
    legend.position  = "right",
    plot.title       = element_text(size = 5, hjust = 0.5)
  ) +
  stat_pwc(
    aes(group = subtype),
    tip.length       = 0,
    size        = 0.116,
    label.size       = 5 / .pt,
    method           = "wilcox.test",
    p.adjust.method  = "BH",
    hide.ns          = TRUE,
    label            = "p.adj.signif",
    y.position       = 0,
    vjust            = 0.5,
    step.increase    = 0.02
  ) +
  labs(x = NULL, y = "G2M score")

ggsave(
  filename = file.path(fig_dir, "figS3c_left.pdf"),
  plot = egg::set_panel_size(
    g2m_score_boxplot,
    width  = grid::unit(1.5, "cm"),
    height = grid::unit(2.0, "cm")
  ),
  width = 6, height = 6, units = "cm"
)

# ---- per-sample cell cycle phase proportions ----
phase_prop_df <- seurat_bmal@meta.data %>%
  group_by(orig.ident, Phase) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  tidyr::complete(orig.ident, Phase, fill = list(count = 0, proportion = 0)) %>%
  left_join(
    seurat_bmal@meta.data %>%
      dplyr::select(orig.ident, subtype) %>%
      distinct(),
    by = "orig.ident"
  )

# ---- plot: phase proportion by subtype (boxplot) ----
g2m_prop_boxplot <- ggplot(phase_prop_df, aes(x = subtype, y = proportion, fill = subtype)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.116, width = 0.6) +
  scale_fill_manual(values = subtype_col1) +
  theme_minimal(base_size = 5) +
  theme(
    panel.grid       = element_blank(),
    panel.background = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.116),
    axis.ticks       = element_line(colour = "black", linewidth = 0.116),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(colour = "black", size = 5),
    legend.position  = "right",
    plot.title       = element_text(size = 5, hjust = 0.5)
  ) +
  stat_pwc(
    aes(group = subtype),
    tip.length       = 0,
    size        = 0.116,
    label.size       = 5 / .pt,
    method           = "wilcox.test",
    p.adjust.method  = "BH",
    hide.ns          = TRUE,
    label            = "p.adj.signif",
    y.position       = 0,
    vjust            = 0.5,
    step.increase    = 0.02
  ) +
  labs(x = NULL, y = "G2M proportion")

ggsave(
  filename = file.path(fig_dir, "figS3c_right.pdf"),
  plot = egg::set_panel_size(
    g2m_prop_boxplot,
    width  = grid::unit(1.5, "cm"),
    height = grid::unit(2.0, "cm")
  ),
  width = 6, height = 6, units = "cm"
)
