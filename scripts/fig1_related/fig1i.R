# ---- packages ----
library(tidyverse)
library(clusterProfiler)
library(edgeR)
library(enrichplot)
library(ggplot2)
library(ggridges)
library(egg)

# ---- paths & sources ----
data_dir  <- "../../data/bulk_rnaseq"
src_dir   <- "../source"
fig_dir   <- "."

bulk_expr_path     <- file.path(data_dir, "bulk_expression_matrix.txt")
bulk_clinical_path <- file.path(data_dir, "bulk_clinical_metadata.txt")

source(file.path(src_dir, "custom_colors.R"))
genesets <- readRDS(file.path(src_dir, "msigdbr_gene_sets.rds"))

term2gene_df <- genesets[, c("pathID", "geneID")] %>% distinct()
id2name_map  <- genesets %>% distinct(pathID, pathID_rename)
id2name_vec  <- setNames(id2name_map$pathID_rename, id2name_map$pathID)

# ---- clinical data (tumor only) ----
clinical_df <- read.table(
  bulk_clinical_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)
clinical_tumor_df <- clinical_df[clinical_df$Group == "Tumor", , drop = FALSE]

# ---- load expression ----
bulk_exp_mat <- read.table(
  bulk_expr_path, header = TRUE, row.names = 1, check.names = FALSE
)
bulk_exp_mat <- bulk_exp_mat[, rownames(clinical_tumor_df), drop = FALSE]
bulk_exp_mat <- as.matrix(log2(bulk_exp_mat + 1))

# ---- sample groups by subtype ----
lme1_samples <- rownames(clinical_tumor_df)[clinical_tumor_df$LME_subtype == "LME_1"]
lme2_samples <- rownames(clinical_tumor_df)[clinical_tumor_df$LME_subtype == "LME_2"]
lme3_samples <- rownames(clinical_tumor_df)[clinical_tumor_df$LME_subtype == "LME_3"]

exp_c1vc3 <- bulk_exp_mat[, c(lme1_samples, lme3_samples), drop = FALSE]
exp_c2vc3 <- bulk_exp_mat[, c(lme2_samples, lme3_samples), drop = FALSE]

group_c1vc3 <- factor(rep(c("LME_1", "LME_3"),
                          times = c(length(lme1_samples), length(lme3_samples))),
                      levels = c("LME_3", "LME_1"))
group_c2vc3 <- factor(rep(c("LME_2", "LME_3"),
                          times = c(length(lme2_samples), length(lme3_samples))),
                      levels = c("LME_3", "LME_2"))

# ---- DE + GSEA helper ----
run_de_gsea <- function(expr_counts, group_fac, term2gene, out_prefix) {
  design <- model.matrix(~ group_fac)
  dge    <- DGEList(counts = expr_counts)
  dge    <- calcNormFactors(dge, method = "TMM")
  v      <- voom(dge, design, plot = FALSE)
  fit    <- lmFit(v, design)
  fit    <- eBayes(fit)
  tt     <- topTable(fit, number = nrow(v$E), sort.by = "none")
  # Create ranked gene list by logFC (names = gene symbol / rowname)
  gene_rank <- tt$logFC
  names(gene_rank) <- rownames(tt)
  gene_rank <- sort(gene_rank, decreasing = TRUE)
  
  gsea_res <- GSEA(gene_rank,
                   TERM2GENE     = term2gene,
                   pvalueCutoff  = 1,
                   pAdjustMethod = "BH")
  
  if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
    gsea_df <- gsea_res@result
  } else {
    gsea_df <- tibble()
  }
  
  write.csv(gsea_df, file.path(fig_dir, paste0(out_prefix, ".result.csv")), row.names = FALSE)
  return(list(top_table = tt, gsea_df = gsea_df))
}

# ---- run comparisons ----
res_c1vc3 <- run_de_gsea(exp_c1vc3, group_c1vc3, term2gene_df, "C1vC3")
res_c2vc3 <- run_de_gsea(exp_c2vc3, group_c2vc3, term2gene_df, "C2vC3")

# ---- rename pathway IDs to readable names ----
rename_description <- function(gsea_df) {
  if (nrow(gsea_df) == 0) return(gsea_df)
  gsea_df$Description <- unname(id2name_vec[ gsea_df$Description ])
  gsea_df
}
gsea_c1vc3 <- rename_description(res_c1vc3$gsea_df)
gsea_c2vc3 <- rename_description(res_c2vc3$gsea_df)

# ---- ridge plot builder ----
prep_gsea_for_ridge <- function(gsea_df, keep_terms, top_table) {
  if (nrow(gsea_df) == 0) return(tibble())
  df <- gsea_df %>%
    filter(Description %in% keep_terms) %>%
    arrange(pvalue) %>%
    mutate(log10P = -log10(pvalue)) %>%
    separate_rows(core_enrichment, sep = "/")
  gene_list <- tibble(gene = rownames(top_table), logfc = top_table$logFC)
  df %>%
    left_join(gene_list, by = c("core_enrichment" = "gene")) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))
}

# ---- C1 vs C3 ridge plot ----
keep_terms_c1vc3 <- c(
  "Angiogenesis [HALLMARK]",
  "Regulation of endothelial cell development [GOBP]",
  "Regulation of endothelial cell differentiation [GOBP]",
  "Tgf beta signaling [HALLMARK]",
  "Extracellular matrix assembly [GOBP]",
  "Extracellular matrix constituent secretion [GOBP]"
)
gsea_c1vc3_sub <- prep_gsea_for_ridge(gsea_c1vc3, keep_terms_c1vc3, res_c1vc3$top_table)

col_scale_acton <- scico::scico(10, direction = -1, palette = "acton")
col_scale_mako  <- viridis::mako(13, direction = -1)[1:10]
custom_colors   <- c(rev(col_scale_mako[1:5]), col_scale_acton[1:8])

p_c1vc3 <- ggplot(gsea_c1vc3_sub, aes(x = logfc, y = Description, fill = log10P)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5, size = 0.25 * 25.4 / 72) +
  labs(x = "log2(fold change)", y = "") +
  scale_fill_gradientn(colors = custom_colors) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25 * 25.4 / 72),
    axis.line = element_line(colour = "black", size = 0.25 * 25.4 / 72),
    axis.ticks = element_line(colour = "black", size = 0.25 * 25.4 / 72),
    axis.text = element_text(color = "black", size = 5),
    axis.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5)
  )

ggsave(
  file.path(fig_dir, "fig1i_left.pdf"),
  egg::set_panel_size(p_c1vc3, width = unit(1.5, "cm"), height = unit(2, "cm")),
  height = 8, width = 16, units = "cm"
)

# ---- C2 vs C3 ridge plot ----
keep_terms_c2vc3 <- c(
  "Antigen processing and presentation [GOBP]",
  "Mhc protein complex [GOCC]",
  "Interferon alpha response [HALLMARK]",
  "T cell receptor signaling pathway [GOBP]",
  "B cell receptor signaling pathway [GOBP]",
  "B cell activation involved in immune response [GOBP]"
)
gsea_c2vc3_sub <- prep_gsea_for_ridge(gsea_c2vc3, keep_terms_c2vc3, res_c2vc3$top_table)

p_c2vc3 <- ggplot(gsea_c2vc3_sub, aes(x = logfc, y = Description, fill = log10P)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5, size = 0.25 * 25.4 / 72) +
  labs(x = "log2(fold change)", y = "") +
  scale_fill_gradientn(colors = custom_colors) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.25 * 25.4 / 72),
    axis.line = element_line(colour = "black", size = 0.25 * 25.4 / 72),
    axis.ticks = element_line(colour = "black", size = 0.25 * 25.4 / 72),
    axis.text = element_text(color = "black", size = 5),
    axis.title = element_text(size = 5),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5)
  )

ggsave(
  file.path(fig_dir, "fig1i_right.pdf"),
  egg::set_panel_size(p_c2vc3, width = unit(1.5, "cm"), height = unit(2, "cm")),
  height = 8, width = 16, units = "cm"
)
