# ---- packages ----
library(Seurat)
library(AUCell)
library(CytoTRACE2)
library(tidyverse)
library(ggplot2)
library(egg)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir <- "../source"
fig_dir <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- inputs ----
B_cell_seurat <- readRDS(file.path(data_dir, "snrna_b_subset_seurat.rds"))
B_Mal_seurat <- subset(B_cell_seurat, fine_cell_type == "B_Mal")

expr_counts <- as.matrix(B_Mal_seurat@assays$RNA@counts)

# ---- AUCell: rankings & AUC ----
cell_rankings <- AUCell_buildRankings(expr_counts, nCores = 12, plotStats = TRUE)

gene_sets_df <- read.csv(file.path(src_dir, "MP_list.csv"),
                         header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gene_sets <- as.list(gene_sets_df)
names(gene_sets)

cells_auc <- AUCell_calcAUC(gene_sets, cell_rankings)
saveRDS(cells_auc, "cells_AUC.rds")

auc_mat <- AUCell::getAUC(cells_auc)  # genesets x cells
stopifnot(identical(rownames(B_Mal_seurat@meta.data), colnames(auc_mat)))

B_Mal_seurat@meta.data <- cbind(B_Mal_seurat@meta.data, t(auc_mat))

# ---- correlation matrix among gene sets (Spearman) ----
auc_mat <- AUCell::getAUC(cells_auc)               # genesets x cells
cells_auc_df <- as.data.frame(t(auc_mat))          # cells x genesets

cor_mat <- cor(cells_auc_df, method = "spearman")

# lower-triangle long-format without duplicates
cor_df_full <- as.data.frame(as.table(cor_mat)) |>
  setNames(c("Var1", "Var2", "Correlation")) |>
  dplyr::filter(Var1 != Var2)

# keep only lower triangle based on factor order
f1 <- as.numeric(factor(cor_df_full$Var1, levels = colnames(cor_mat)))
f2 <- as.numeric(factor(cor_df_full$Var2, levels = colnames(cor_mat)))
cor_df <- cor_df_full[f1 < f2, , drop = FALSE]

# ---- plot: lower-triangle bubble plot ----
max_abs <- max(abs(cor_df$Correlation), na.rm = TRUE)

cor_dotplot <- ggplot(cor_df, aes(x = Var1, y = Var2, size = abs(Correlation), color = Correlation)) +
  geom_point() +
  scale_y_discrete(limits = rev) +
  scale_color_gradient2(low = "#1d91c0", mid = "white", high = "#91003f",
                        midpoint = 0, limits = c(-max_abs, max_abs)) +
  scale_size(range = c(1, 2), limits = c(0, max_abs)) +
  theme_bw(base_size = 5) +
  labs(x = NULL, y = NULL) +
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
    fill = guide_colourbar(title.position = "top", title.hjust = 0.5)
  )

ggsave(
  filename = file.path(fig_dir, "fig3c.pdf"),
  plot = egg::set_panel_size(
    cor_dotplot,
    width  = grid::unit(2, "cm"),
    height = grid::unit(2, "cm")
  ),
  width = 10, height = 10, units = "cm"
)

# ---- CytoTRACE2 integration & correlation with AUCell scores ----
cytotrace_sce <- cytotrace2(B_Mal_seurat, 
                            is_seurat = TRUE, 
                            slot_type = "counts", 
                            species = 'human',
                            seed = 421)
stopifnot(identical(rownames(B_Mal_seurat@meta.data), rownames(cytotrace_sce@meta.data)))

B_Mal_seurat@meta.data <- cbind(
  B_Mal_seurat@meta.data,
  cytotrace_sce@meta.data[, c("CytoTRACE2_Score", "CytoTRACE2_Potency")]
)

gene_set_names <- names(gene_sets)
meta_df <- B_Mal_seurat@meta.data

cor_results <- sapply(gene_set_names, function(col) {
  test <- cor.test(meta_df[["CytoTRACE2_Score"]], meta_df[[col]], method = "spearman")
  c(correlation = unname(test$estimate), p_value = unname(test$p.value))
})

lollipop_df <- data.frame(
  MP  = colnames(cor_results),
  rho = cor_results["correlation", ]
) |>
  arrange(rho) |>
  mutate(MP = factor(MP, levels = MP))

lollipop_plot <- ggplot(lollipop_df, aes(x = rho, y = MP)) +
  geom_segment(aes(y = MP, yend = MP, x = 0, xend = rho), color = "black", linewidth = 0.232) +
  geom_point(size = 1, pch = 19, color = "#7f9db4") +
  xlab("Correlation coefficient") + ylab(NULL) +
  theme_bw(base_size = 5) +
  theme(
    text = element_text(size = 5),
    panel.grid.major = element_line(colour = "grey90", size = 0.116),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", size = 0.116),
    axis.line = element_line(colour = "black", size = 0.116),
    axis.ticks = element_line(colour = "black", size = 0.116),
    axis.text.x = element_text(colour = "black", size = 5),
    axis.text.y = element_text(colour = "black", size = 5)
  )

ggsave(
  filename = file.path(fig_dir, "fig3d.pdf"),
  plot = egg::set_panel_size(
    lollipop_plot,
    width  = grid::unit(1.5, "cm"),
    height = grid::unit(2.0, "cm")
  ),
  width = 6, height = 6, units = "cm"
)