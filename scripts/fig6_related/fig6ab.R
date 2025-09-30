# ---- packages ----
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(egg)
library(grid)
library(ComplexHeatmap)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- inputs ----
seurat_all <- readRDS(file.path(data_dir, "snrna_all_cells_seurat.rds"))
seurat_use <- subset(
  seurat_all,
  fine_cell_type %in% c("Fibro_TCF21", "AT1_AGER", "AT2_SFTPC", 
                        "Hepa_ALB", "Endo_HPGD", "B_Mal", "B_Nor"),
  invert = TRUE
)
seurat_use$coarse_cell_type <- droplevels(seurat_use$coarse_cell_type)
seurat_use$fine_cell_type   <- droplevels(seurat_use$fine_cell_type)

meta_df <- seurat_use@meta.data

# ---- compute per-coarse cell type proportions ----
coarse_types <- unique(meta_df$coarse_cell_type)
prop_results <- NULL

for (coarse_ct in coarse_types) {
  
  subset_meta <- meta_df[meta_df$coarse_cell_type == coarse_ct, ]
  
  subset_meta$orig.ident     <- factor(subset_meta$orig.ident, 
                                       levels = unique(meta_df$orig.ident))
  subset_meta$fine_cell_type <- droplevels(subset_meta$fine_cell_type)
  
  prop_mat <- prop.table(table(subset_meta$orig.ident, 
                               subset_meta$fine_cell_type), margin = 1)
  
  if (is.null(prop_results)) {
    prop_results <- prop_mat
  } else {
    prop_results <- cbind(prop_results, prop_mat)
  }
}

# ---- correlation analysis ----
cor_mat <- cor(prop_results, method = "spearman")

cor_test_mat <- function(df) {
  n <- ncol(df)
  p_mat <- matrix(NA, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(df[, i], df[, j], method = "spearman")
      p_mat[i, j] <- test$p.value
      p_mat[j, i] <- test$p.value
    }
  }
  diag(p_mat) <- 0
  return(p_mat)
}

p_mat <- cor_test_mat(prop_results)

# ---- heatmap ----
pdf('fig6a.pdf', width = 8, height = 8)
Heatmap(cor_mat,
        name = "Spearman\ncorrelation",
        col = heatmap_col5,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_dend_gp = gpar(lwd = 0.33),
        row_dend_reorder = FALSE,
        column_dend_reorder = FALSE,
        show_row_dend = TRUE,
        show_column_dend = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        width = unit(ncol(cor_mat) * 0.24, 'cm'),
        height = unit(nrow(cor_mat) * 0.24, 'cm'),
        column_title_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5),
        row_title_gp = gpar(fontsize = 5),
        row_names_gp = gpar(fontsize = 5))
dev.off()

# ---- grouped proportions (CM1 vs CM2) ----
df_prop <- prop.table(table(meta_df$orig.ident, meta_df$fine_cell_type), margin = 1)

df_sum <- data.frame(
  CM1 = rowSums(df_prop[, c("cDC_XCR1", "pDC_CLEC4C", "FDC_CR2", "Mac_SPARC", 
                            "Fibro_CCL19", "CD4T_CXCL13", "CD8T_HAVCR2", 
                            "CD4T_FOXP3", "cDC_CLEC10A", "Endo_CCL21")]),
  CM2 = rowSums(df_prop[, c("Mac_MARCO", "Endo_ACKR1", "CD8T_GNLY", "Mac_STAB1", 
                            "CD8T_CCL5", "Neu_G0S2", "Endo_CXCL12", "Peri_PDGFRB")])
)

df_sum$subtype <- sapply(
  rownames(df_sum),
  function(x) unique(meta_df[meta_df$orig.ident == x, "subtype"])
)

df_melt <- reshape2::melt(df_sum, id.vars = "subtype")
df_melt$subtype <- factor(df_melt$subtype, levels = c("LME_1", "LME_2", "LME_3"))

plot_prop <- ggplot(df_melt, aes(x = subtype, y = value)) +
  geom_jitter(aes(fill = subtype), width = 0.2, alpha = 0.8, size = 0.75,
              stroke = 0.116, shape = 21) +
  geom_boxplot(outlier.shape = NA, size = 0.116, width = 0.3, 
               color = "black", fill = NA) +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = "", y = "Proportion (%)") +
  scale_fill_manual(values = subtype_col1) +
  theme(text = element_text(size = 5),
        strip.text = element_text(size = 5),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.116),
        axis.ticks = element_line(colour = "black", size = 0.116),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 5),
        plot.title = element_text(size = 5, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 5)) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.5),
    labels = function(x) x * 100
  ) +
  ggpubr::stat_pwc(aes(group = subtype),
                   label.size = 5 * 25.4 / 72,
                   tip.length = 0,
                   size = 0.116,
                   hide.ns = TRUE,
                   vjust = 0.5,
                   method = "wilcox.test",
                   p.adjust.method = "BH",
                   label = "p.adj.signif") +
  coord_cartesian(clip = "off")

ggsave(filename = 'fig6b.pdf', width = 8, height = 5, units = 'cm',
       egg::set_panel_size(plot_prop, width = unit(1.5, "cm"),
                           height = unit(1.5, "cm")))