# ---- packages ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(egg)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))  # uses subtype_col2

# ---- all cells ----
all_cells_seurat <- readRDS(file.path(data_dir, "snrna_all_cells_seurat.rds"))

use_cells_seurat <- subset(
  all_cells_seurat,
  fine_cell_type %in% c("Fibro_TCF21", "AT1_AGER", "AT2_SFTPC", "Hepa_ALB", "Endo_HPGD"),
  invert = TRUE
)
use_cells_seurat$fine_cell_type <- droplevels(use_cells_seurat$fine_cell_type)

meta_df <- use_cells_seurat@meta.data

# ---- odds ratio & p-value matrices (cell type × subtype) ----
fine_types <- sort(unique(meta_df$fine_cell_type))
subtypes   <- sort(unique(meta_df$subtype))

or_matrix <- matrix(NA_real_,
                    nrow = length(fine_types),
                    ncol = length(subtypes),
                    dimnames = list(fine_types, subtypes))

p_matrix  <- or_matrix

for (cell_i in rownames(or_matrix)) {
  for (subtype_j in colnames(or_matrix)) {
    tbl <- table(meta_df$fine_cell_type == cell_i, meta_df$subtype == subtype_j)
    ft  <- fisher.test(tbl)
    or_matrix[cell_i, subtype_j] <- unname(ft$estimate)
    p_matrix[cell_i, subtype_j]  <- ft$p.value
  }
}

padj_matrix <- matrix(
  p.adjust(as.vector(p_matrix), method = "BH"),
  nrow = nrow(p_matrix),
  ncol = ncol(p_matrix),
  dimnames = dimnames(p_matrix)
)

merged_matrix <- cbind(round(padj_matrix, 4), round(or_matrix, 2))

# ---- long-format for plotting ----
# OR long
or_df <- or_matrix %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("fine_cell_type") %>%
  pivot_longer(-fine_cell_type, names_to = "subtype", values_to = "OR")

# adjusted p long
padj_df <- padj_matrix %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("fine_cell_type") %>%
  pivot_longer(-fine_cell_type, names_to = "subtype", values_to = "p_adjusted")

# annotate coarse cell type by pattern
or_df <- or_df %>%
  mutate(
    coarse_cell_type = case_when(
      grepl("^B_",        fine_cell_type) ~ "B_cell",
      grepl("^(CD8T_|CD4T_)", fine_cell_type) ~ "T_cell",
      grepl("^(Mac_|cDC_|pDC_|Neu_|Mast_)", fine_cell_type) ~ "Myeloid_cell",
      grepl("^(Endo_|Peri_)",     fine_cell_type) ~ "Endo_cell",
      grepl("^(Fibro_|FDC_)", fine_cell_type) ~ "Fibro",
      TRUE ~ "Other"
    )
  )

plot_df <- or_df %>%
  left_join(padj_df, by = c("fine_cell_type", "subtype")) %>%
  mutate(
    sig = if_else(p_adjusted < 0.05, "sig", "no.sig")
  )
plot_df$fine_cell_type <- factor(plot_df$fine_cell_type, levels = levels(fine_types))

# ---- plot ----
or_plot <- ggplot(
  plot_df,
  aes(
    x = fine_cell_type,
    y = log(OR),
    shape = sig,
    color = subtype,
    group = interaction(subtype, coarse_cell_type)
  )
) +
  geom_line(size = 0.116) +
  geom_point(size = 1, stroke = 0.116) +
  scale_shape_manual(values = c("no.sig" = 21, "sig" = 16)) +
  scale_color_manual(values = subtype_col2) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.116) +
  labs(x = "", y = "ln(odds ratio)") +
  theme(
    text = element_text(size = 5),
    panel.grid.major.x = element_line(colour = "grey90", size = 0.116),
    panel.grid.minor   = element_blank(),
    panel.background   = element_blank(),
    axis.line          = element_line(colour = "black", size = 0.116),
    axis.ticks.x       = element_line(colour = "black", size = 0.116),
    axis.ticks.y       = element_line(colour = "black", size = 0.116),
    axis.title         = element_text(colour = "black", size = 5),
    axis.text.x        = element_text(colour = "black", angle = 60, vjust = 1, hjust = 1, size = 5),
    axis.text.y        = element_text(colour = "black", size = 5),
    legend.position    = "bottom",
    legend.title       = element_text(size = 5),
    legend.text        = element_text(size = 5)
  )

ggsave(
  file.path(fig_dir, "fig2f.pdf"),
  egg::set_panel_size(or_plot, width = unit(6, "cm"), height = unit(2, "cm")),
  width = 12, height = 6, units = "cm"
)