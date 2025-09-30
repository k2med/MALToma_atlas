# ---- packages ----
library(edgeR)
library(limma)
library(dplyr)
library(ggvenn)

# ---- paths & sources ----
data_dir  <- "../../data/bulk_rnaseq"
src_dir   <- "../source"
fig_dir   <- "."

bulk_expr_path     <- file.path(data_dir, "bulk_expression_matrix.txt")
bulk_clinical_path <- file.path(data_dir, "bulk_clinical_metadata.txt")

source(file.path(src_dir, "custom_colors.R")) 

# ---- data ----
bulk_expr_mat <- read.table(bulk_expr_path, header = TRUE, row.names = 1, check.names = FALSE)
bulk_expr_mat <- as.matrix(log2(bulk_expr_mat + 1))

clinical_df <- read.table(
  bulk_clinical_path, header = TRUE, sep = "\t", row.names = 1, na.strings = c("", "NA")
)

lme1_samples   <- rownames(clinical_df)[which(clinical_df$LME_subtype == "LME_1")]
lme2_samples   <- rownames(clinical_df)[which(clinical_df$LME_subtype == "LME_2")]
lme3_samples   <- rownames(clinical_df)[which(clinical_df$LME_subtype == "LME_3")]
normal_samples <- rownames(clinical_df)[which(clinical_df$Group == "Normal")]

bulk_expr_mat_LME1vNormal <- bulk_expr_mat[, c(normal_samples, lme1_samples), drop = FALSE]
bulk_expr_mat_LME2vNormal <- bulk_expr_mat[, c(normal_samples, lme2_samples), drop = FALSE]
bulk_expr_mat_LME3vNormal <- bulk_expr_mat[, c(normal_samples, lme3_samples), drop = FALSE]

group_LME1vNormal <- factor(
  rep(c("Normal", "LME_1"), times = c(length(normal_samples), length(lme1_samples))),
  levels = c("Normal", "LME_1")
)
group_LME2vNormal <- factor(
  rep(c("Normal", "LME_2"), times = c(length(normal_samples), length(lme2_samples))),
  levels = c("Normal", "LME_2")
)
group_LME3vNormal <- factor(
  rep(c("Normal", "LME_3"), times = c(length(normal_samples), length(lme3_samples))),
  levels = c("Normal", "LME_3")
)

# ---- DE loops ----
for (compare_i in c("LME1vNormal", "LME2vNormal", "LME3vNormal")) {
  
  counts_tmp <- get(paste0("bulk_expr_mat_", compare_i))
  group_tmp  <- get(paste0("group_", compare_i))
  
  design <- model.matrix(~ group_tmp)
  
  dge <- DGEList(counts = counts_tmp)
  dge <- calcNormFactors(dge, method = "TMM")
  
  v   <- voom(dge, design, plot = TRUE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  gene_diff  <- topTable(fit, number = nrow(v$E))
  
  gene_diff[gene_diff$logFC >=  2 & gene_diff$adj.P.Val < 0.05, "sig"] <- "up"
  gene_diff[gene_diff$logFC <= -2 & gene_diff$adj.P.Val < 0.05, "sig"] <- "down"
  gene_diff[abs(gene_diff$logFC) <= 2 | gene_diff$adj.P.Val >= 0.05, "sig"] <- "none"
  
  gene_diff_up   <- subset(gene_diff, sig == "up")
  
  write.table(gene_diff_up, paste0(compare_i, ".up.txt"),
              sep = "\t", col.names = NA, quote = FALSE)
}

# ---- Venn ----
LME1_up <- read.table("LME1vNormal.up.txt")
LME2_up <- read.table("LME2vNormal.up.txt")
LME3_up <- read.table("LME3vNormal.up.txt")

gene_list <- list(
  LME1 = rownames(LME1_up),
  LME2 = rownames(LME2_up),
  LME3 = rownames(LME3_up)
)

ggvenn(
  gene_list,
  show_percentage = FALSE,
  stroke_color    = "white",
  fill_color      = unname(subtype_col1),
  set_name_color  = unname(subtype_col2),
  stroke_size     = 0.25,
  text_size       = 5 * 25.4 / 72,
  set_name_size   = 5 * 25.4 / 72
)
ggsave("figS1f.pdf", height = 3, width = 3, units = "cm")