# ---- packages ----
library(ggplot2)
library(dplyr)
library(egg)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/bulk_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

bulk_expr_path     <- file.path(data_dir, "bulk_expression_matrix.txt")
bulk_clinical_path <- file.path(data_dir, "bulk_clinical_metadata.txt")

source(file.path(src_dir, "custom_colors.R"))

# ---- load expression & prepare PCA input ----
bulk_expr_mat <- read.table(
  bulk_expr_path, header = TRUE, row.names = 1, check.names = FALSE
)

pca_input_df <- t(log2(bulk_expr_mat + 1))
pca_input_df <- pca_input_df[, apply(pca_input_df, 2, var) != 0, drop = FALSE]

# ---- PCA fit & scores ----
pca_fit <- prcomp(pca_input_df, center = TRUE, scale. = TRUE)
pca_scores_mat <- predict(pca_fit)
pca_scores_df  <- data.frame(
  PC1 = pca_scores_mat[, 1],
  PC2 = pca_scores_mat[, 2],
  check.names = FALSE
)

pca_var_exp_vec <- summary(pca_fit)$importance["Proportion of Variance", ]

# ---- load clinical metadata ----
clinical_df <- read.table(
  bulk_clinical_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)

# standard clinical variables to keep
clinical_keep_vec <- c(
  "Cohort", "Sex", "Age", "Group", "Site", "Ann_Arbor_stage", "LME_subtype", "BIRC3_MALT_fusion"
)

# merge with PCA scores
pca_scores_df <- cbind(
  pca_scores_df,
  clinical_df[rownames(pca_scores_df),
              intersect(clinical_keep_vec, colnames(clinical_df)),
              drop = FALSE]
)

# ---- plotting function ----
plot_pca <- function(data, group_col = "cohort", color_pal,
                        var_exp = c(0.5, 0.3),
                        out_file = file.path(fig_dir, paste0("pca_by_", group_col, ".pdf"))) {
  
  p <- ggplot(data = data, aes(x = PC1, y = PC2)) +
    geom_point(aes_string(color = group_col), size = 0.1) +
    scale_color_manual(values = color_pal) +
    labs(
      x = paste0("PC1 (", round(var_exp[1], 3) * 100, "%)"),
      y = paste0("PC2 (", round(var_exp[2], 3) * 100, "%)")
    ) +
    theme_classic() +
    theme(
      axis.line.x   = element_line(colour = "black", size = 0.116),
      axis.line.y   = element_line(colour = "black", size = 0.116),
      axis.ticks.x  = element_line(colour = "black", size = 0.116),
      axis.ticks.y  = element_line(colour = "black", size = 0.116),
      axis.text.x   = element_text(color = "black", size = 5),
      axis.text.y   = element_text(color = "black", size = 5),
      axis.title    = element_text(color = "black", size = 5),
      strip.text    = element_blank(),
      legend.text   = element_text(size = 5),
      legend.title  = element_text(size = 5),
      aspect.ratio  = 1
    )
  
  ggsave(
    out_file,
    plot  = egg::set_panel_size(p, width = unit(2, "cm"), height = unit(2, "cm")),
    width = 10, height = 10, units = "cm"
  )
}

# ---- run plots ----
plot_pca(
  data      = pca_scores_df,
  group_col = "Cohort",
  color_pal = cohort_col,
  var_exp   = pca_var_exp_vec,
  out_file  = "fig1b.pdf"
)

plot_pca(
  data      = pca_scores_df,
  group_col = "Group",
  color_pal = group_col,
  var_exp   = pca_var_exp_vec,
  out_file  = "fig1c.pdf"
)

plot_pca(
  data      = pca_scores_df %>%
    mutate(Site = ifelse(Site %in% c("Ocular_adnexa",
                                     "Stomach",
                                     "Lung",
                                     "Intestine",
                                     "Nasopharynx",
                                     "Lymph_node"),
                         Site, "Others")),
  group_col = "Site",
  color_pal = site_col,
  var_exp   = pca_var_exp_vec,
  out_file  = "figS1b_left.pdf"
)

plot_pca(
  data      = pca_scores_df,
  group_col = "LME_subtype",
  color_pal = subtype_col1,
  var_exp   = pca_var_exp_vec,
  out_file  = "figS1b_right.pdf"
)