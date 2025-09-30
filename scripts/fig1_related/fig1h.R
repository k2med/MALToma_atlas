# ---- packages ----
library(ggtern)
library(ggrastr)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ---- paths & sources ----
data_dir  <- "../../data/bulk_rnaseq"
src_dir   <- "../source"
fig_dir   <- "."

bulk_expr_path     <- file.path(data_dir, "bulk_expression_matrix.txt")
bulk_clinical_path <- file.path(data_dir, "bulk_clinical_metadata.txt")

source(file.path(src_dir, "custom_colors.R"))

# ---- clinical data (tumor only) ----
clinical_df <- read.table(
  bulk_clinical_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)
clinical_tumor_df <- clinical_df[clinical_df$Group == "Tumor", , drop = FALSE]

# ---- load expression ----
bulk_expr_mat <- read.table(
  bulk_expr_path, header = TRUE, row.names = 1, check.names = FALSE
)
bulk_expr_mat <- bulk_expr_mat[, rownames(clinical_tumor_df), drop = FALSE]

# ---- per-subtype averages ----
idx_lme1 <- rownames(clinical_tumor_df)[clinical_tumor_df$LME_subtype == "LME_1"]
idx_lme2 <- rownames(clinical_tumor_df)[clinical_tumor_df$LME_subtype == "LME_2"]
idx_lme3 <- rownames(clinical_tumor_df)[clinical_tumor_df$LME_subtype == "LME_3"]

avg_exp_df <- data.frame(
  LME_1  = rowMeans(as.matrix(bulk_expr_mat[, idx_lme1, drop = FALSE]), na.rm = TRUE),
  LME_2  = rowMeans(as.matrix(bulk_expr_mat[, idx_lme2, drop = FALSE]), na.rm = TRUE),
  LME_3  = rowMeans(as.matrix(bulk_expr_mat[, idx_lme3, drop = FALSE]), na.rm = TRUE),
  gene   = rownames(bulk_expr_mat),
  check.names = FALSE
) %>%
  mutate(average = rowMeans(across(c(LME_1, LME_2, LME_3)), na.rm = TRUE))

# ---- highlight gene sets ----
genes_lme1 <- c("COL3A1","COL6A2","COL6A1","COL1A1","VEGFA","CD151","POSTN","SDC1")
genes_lme2 <- c("PTPRC","CXCL13","TRAC","PDCD1","GZMK","B2M","CD74","CD68")
genes_lme3 <- c("MKI67","CDK1","HMGB2","EZH2","ZWINT","E2F1","MCM2","AURKA")

avg_exp_df <- avg_exp_df %>%
  mutate(
    gene_type = case_when(
      gene %in% genes_lme1 ~ "LME_1",
      gene %in% genes_lme2 ~ "LME_2",
      gene %in% genes_lme3 ~ "LME_3",
      TRUE ~ "Other"
    )
  )

# ---- ternary plot ----
p_tern <- ggtern(data = avg_exp_df, aes(x = LME_1, y = LME_2, z = LME_3)) +
  geom_point_rast(aes(size = average), color = "grey80", alpha = 1, show.legend = FALSE) +
  scale_size(range = c(0, 1)) +
  geom_point(
    data  = subset(avg_exp_df, gene_type == "LME_1"),
    size  = 1, shape = 19, color = subtype_col1[1], show.legend = FALSE
  ) +
  geom_point(
    data  = subset(avg_exp_df, gene_type == "LME_2"),
    size  = 1, shape = 19, color = subtype_col1[2], show.legend = FALSE
  ) +
  geom_point(
    data  = subset(avg_exp_df, gene_type == "LME_3"),
    size  = 1, shape = 19, color = subtype_col1[3], show.legend = FALSE
  ) +
  geom_text(
    data  = subset(avg_exp_df, gene_type != "Other"),
    aes(label = gene),
    color = "black", size = 3 / .pt, fontface = "italic", show.legend = FALSE
  ) +
  theme_bw() +
  theme(
    axis.text  = element_text(color = "black", size = 5),
    axis.title = element_text(size = 5)
  )

ggsave(
  filename = file.path(fig_dir, "fig1h.pdf"),
  plot = p_tern,
  height = 5, width = 5, units = "cm"
)
