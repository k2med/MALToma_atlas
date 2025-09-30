# ---- packages ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(viridis)
library(scico)
library(grid)
library(egg)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- input: factor matrix (cell-type × spatial factors) ----
factor_fp <- file.path(
  "./colocation/cell2location_map/CoLocatedComb",
  "CoLocatedGroupsSklearnNMF_30515locations_25factors",
  "cell_type_fractions_mean",
  "n_fact5.csv"
)
factor_df <- read.csv(factor_fp, header = TRUE, row.names = 1, check.names = FALSE)
colnames(factor_df) <- paste0("SF", 1:5)

# ---- clustering order for rows/cols (for consistent dotplot axes) ----
factor_mat <- as.matrix(factor_df)
pheat <- pheatmap(factor_mat, silent = TRUE)
clustered_cell_types <- rownames(factor_mat)[pheat$tree_row$order]
clustered_factors <- colnames(factor_mat)[pheat$tree_col$order]

# ---- long-format data with CM group annotation ----
factor_long <- factor_df %>%
  tibble::rownames_to_column("cell_type") %>%
  pivot_longer(
    cols = starts_with("SF"),
    names_to = "SF",
    values_to = "Importance"
  )

factor_long$cell_type <- factor(factor_long$cell_type, levels = rev(clustered_cell_types))

# ---- dot plot (cell types × spatial factors) ----
factor_dotplot <- ggplot(factor_long, aes(x = SF, y = cell_type, size = Importance, color = Importance)) +
  geom_point() +
  scale_color_gradientn(colors = scico(10, direction = -1, palette = "acton")[1:8]) +
  scale_size(range = c(0, 2), limits = c(0, 1)) +
  labs(x = NULL, y = NULL, size = "Importance", color = "Importance") +
  theme_bw(base_size = 5) +
  theme(
    text = element_text(size = 5),
    panel.grid.major = element_line(colour = "grey90", size = 0.116),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", size = 0.116),
    axis.line = element_line(colour = "black", size = 0.116),
    axis.ticks = element_line(colour = "black", size = 0.116),
    axis.text.x = element_text(colour = "black", angle = 60, vjust = 1, hjust = 1, size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    legend.position = "right",
    legend.title = element_text(size = 5)
  ) +
  guides(
    size = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(stroke = 0.4)),
    color = guide_colourbar(title.position = "top", title.hjust = 0.5)
  )

ggsave(
  filename = file.path(fig_dir, "fig7b.pdf"),
  plot = egg::set_panel_size(
    factor_dotplot,
    width  = unit(length(unique(factor_long$SF)) * 0.3, "cm"),
    height = unit(length(unique(factor_long$cell_type)) * 0.2, "cm")
  ),
  width = 10, height = 10, units = "cm"
)

# ---- input: per-spot factor loadings ----
spot_fp <- file.path(
  "./colocation_250711//cell2location_map/CoLocatedComb",
  "CoLocatedGroupsSklearnNMF_30515locations_25factors",
  "location_factors_mean",
  "n_fact5.csv"
)
spot_factor <- read.csv(spot_fp, header = TRUE, row.names = 1, check.names = FALSE)
colnames(spot_factor) <- paste0("SF", 1:5)

# derive sample and subtype
spot_factor$sample <- sapply(strsplit(rownames(spot_factor), "_"), `[`, 1)
spot_factor <- spot_factor %>%
  mutate(
    subtype = case_when(
      sample %in% c("S31", "S32", "S71") ~ "LME_1",
      sample %in% c("S01", "S02", "S89") ~ "LME_2",
      sample %in% c("S04", "S06", "S87") ~ "LME_3",
      TRUE ~ NA_character_
    )
  )

# mean factor importance by subtype
spot_factor_mean <- spot_factor %>%
  group_by(subtype) %>%
  summarise(across(SF1:SF5, ~ mean(.x, na.rm = TRUE)), .groups = "drop")

plot_data <- spot_factor_mean %>%
  pivot_longer(
    cols = starts_with("SF"),
    names_to = "SF",
    values_to = "Value"
  ) %>%
  rename(LME = subtype)

# ---- line plot across factors (per LME) ----
line_plot <- ggplot(plot_data, aes(x = SF, y = Value, color = LME, group = LME)) +
  geom_line(size = 0.116) +
  geom_point(size = 1, stroke = 0.116) +
  scale_color_manual(values = subtype_col2) +
  scale_y_log10() +
  labs(x = NULL, y = "Mean importance", title = NULL) +
  theme_minimal(base_size = 5) +
  theme(
    legend.title = element_blank(),
    text = element_text(size = 5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", size = 0.116),
    axis.ticks = element_line(colour = "black", size = 0.116),
    axis.text.x = element_text(colour = "black", size = 5, angle = 60, vjust = 1, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 5),
    legend.position = "right"
  ) +
  coord_cartesian(clip = "off")

ggsave(
  filename = file.path(fig_dir, "fig7d.pdf"),
  plot = egg::set_panel_size(line_plot, width = unit(2.4, "cm"), height = unit(1.2, "cm")),
  width = 10, height = 10, units = "cm"
)