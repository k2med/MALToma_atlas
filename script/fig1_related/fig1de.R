# ---- packages ----
library(tidyverse)
library(clusterProfiler)
library(edgeR)
library(limma)
library(ggrepel)
library(DOSE)
library(org.Hs.eg.db)
library(egg)

# ---- paths & sources ----
data_dir  <- "../../data/bulk_rnaseq"
src_dir   <- "../source"
fig_dir   <- "."

bulk_expr_path     <- file.path(data_dir, "bulk_expression_matrix.txt")
bulk_clinical_path <- file.path(data_dir, "bulk_clinical_metadata.txt")
genesets_path      <- file.path(src_dir, "msigdbr_gene_sets.rds")

source(file.path(src_dir, "custom_colors.R"))
genesets <- readRDS(genesets_path)

# ---- load data ----
clinical_df <- read.table(bulk_clinical_path, header = TRUE, sep = "\t", 
                          row.names = 1, check.names = FALSE)
bulk_exp_mat <- read.table(bulk_expr_path, header = TRUE, row.names = 1, 
                           check.names = FALSE)
bulk_exp_mat <- log2(bulk_exp_mat + 1)

# ---- groups ----
normal_samples <- rownames(clinical_df)[clinical_df$Group == "Normal"]
tumor_samples  <- rownames(clinical_df)[clinical_df$Group == "Tumor"]
group_factor   <- factor(
  rep(c("Normal", "Tumor"), times = c(length(normal_samples), length(tumor_samples))),
  levels = c("Normal", "Tumor")
)
bulk_exp_mat <- bulk_exp_mat[, c(normal_samples, tumor_samples)]

# ---- differential expression ----
design <- model.matrix(~ group_factor)
dge    <- DGEList(counts = bulk_exp_mat)
dge    <- calcNormFactors(dge, method = "TMM")
v      <- voom(dge, design, plot = FALSE)
fit    <- eBayes(lmFit(v, design))
gene_diff <- topTable(fit, number = nrow(v$E))

# ---- helper plotting function ----
plot_ranked <- function(df, value_col, highlight, xlab, ylab, fontface, filename,
                        point.padding, box.padding, max.overlaps) {
  df <- df[order(df[[value_col]], decreasing = TRUE), ]
  df$Index <- seq_len(nrow(df))
  df$Highlight <- ifelse(rownames(df) %in% highlight, "Highlight", "Other")
  
  p <- ggplot(df, aes(x = Index, y = .data[[value_col]])) +
    geom_line(color = "black", size = 0.116) +
    geom_point(
      data = subset(df, Highlight == "Highlight"),
      size = 1, shape = 21, fill = "#f7931e", color = "black", stroke = 0.116
    ) +
    geom_text_repel(
      data = subset(df, Highlight == "Highlight"),
      aes(label = rownames(df)[Index]),
      color = "black", size = 5 / .pt, segment.size = 0.116,
      min.segment.length = 0.1, point.padding = point.padding, box.padding = box.padding,
      max.overlaps = max.overlaps, fontface = fontface
    ) +
    theme_minimal(base_size = 5) +
    theme(
      panel.grid = element_blank(),
      axis.line  = element_line(colour = "black", size = 0.116),
      axis.ticks = element_line(colour = "black", size = 0.116),
      axis.text  = element_text(colour = "black", size = 5),
      legend.position = "none"
    ) +
    labs(x = xlab, y = ylab, title = "") +
    coord_cartesian(clip = "off")
  
  ggsave(file.path(fig_dir, filename),
         egg::set_panel_size(p, width = unit(3.2, "cm"), height = unit(1.6, "cm")),
         width = 6, height = 6, units = "cm")
}

# ---- plot 1: upregulated genes ----
highlight_genes <- c("CD79A", "MS4A1", "PAX5", "SPN", "CD40")
upregulated <- subset(gene_diff, logFC > 0)
plot_ranked(upregulated, "logFC", highlight_genes,
            xlab = "Upregulated genes ranked by log2FC",
            ylab = "log2FC",
            fontface = "italic",
            point.padding = 0, 
            box.padding = 1, 
            max.overlaps = 20,
            filename = "fig1d.pdf")

# ---- GSEA (HALLMARK) ----
genesets_hallmark <- genesets[grep("HALLMARK", genesets$pathID), ]
gene_list <- gene_diff$logFC
names(gene_list) <- rownames(gene_diff)
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_run <- GSEA(gene_list, TERM2GENE = genesets_hallmark,
                 pvalueCutoff = 1, pAdjustMethod = "BH")
gsea_res <- gsea_run@result[order(gsea_run@result$NES, decreasing = TRUE), ]
gsea_res$Index <- seq_len(nrow(gsea_res))

# ---- plot 2: GSEA results ----
highlight_pathways <- c(
  "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS",
  "HALLMARK_DNA_REPAIR", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
)
rownames(gsea_res) <- gsea_res$ID
plot_ranked(gsea_res, "NES", highlight_pathways,
            xlab = "Pathways ranked by NES",
            ylab = "NES",
            fontface = "plain",
            point.padding = 0, 
            box.padding = 0.3, 
            max.overlaps = 10,
            filename = "fig1e.pdf")