# ---- packages ----
library(Seurat)
library(AnnoProbe)
library(rjags)
library(infercnv)
library(gplots)
library(ggplot2)
library(dplyr)
library(Nebulosa)
library(RColorBrewer)
library(ggridges)
library(ggpubr)
library(AUCell)
library(egg)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- load B-cell Seurat ----
B_cell_seurat <- readRDS(file.path(data_dir, "snrna_b_subset_seurat.rds"))

# ---- inferCNV input files (group & gene order) ----
group_info  <- data.frame(v1 = colnames(B_cell_seurat),
                          v2 = B_cell_seurat$seurat_clusters,
                          check.names = FALSE)
group_file  <- "groupFiles.txt"
write.table(group_info, file = group_file, sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)

gene_info <- annoGene(rownames(B_cell_seurat), "SYMBOL", "human")
gene_info <- gene_info[with(gene_info, order(chr, start)), c(1,4:6)]
gene_info <- gene_info[!duplicated(gene_info[, 1]), ]
gene_info$chr <- factor(
  gene_info$chr,
  levels = c(
    paste0("chr", 1:22), "chrM", "chrX", "chrY"
  )
)
gene_info <- gene_info[order(gene_info$chr), ]
gene_file <- "geneFile.txt"
write.table(gene_info, file = gene_file, sep = "\t", quote = FALSE,
            col.names = FALSE, row.names = FALSE)

# ---- counts matrix aligned to gene order ----
exp_mat  <- as.data.frame(GetAssayData(B_cell_seurat, slot = "counts", assay = "RNA"))
exp_keep <- exp_mat[rownames(exp_mat) %in% gene_info[, 1], , drop = FALSE]
exp_keep <- exp_keep[match(gene_info[, 1], rownames(exp_keep)), , drop = FALSE]

# ---- run inferCNV ----
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix      = exp_keep,
  annotations_file       = group_file,
  delim                  = "\t",
  gene_order_file        = gene_file,
  min_max_counts_per_cell = c(1, +Inf),
  ref_group_names        = NULL
)

infercnv_obj <- run(
  infercnv_obj,
  cutoff            = 0.1,
  out_dir           = "B_output",
  cluster_by_groups = FALSE,
  plot_steps        = FALSE,
  HMM               = FALSE,
  denoise           = TRUE
)

# ---- CNV score per cell ----
cnv_mat <- infercnv_obj@expr.data
cnv_score_df <- data.frame(
  cell = colnames(cnv_mat),
  CNV_score = colMeans((cnv_mat - 1)^2),
  row.names = colnames(cnv_mat),
  check.names = FALSE
)

# join by cell names
B_cell_seurat$CNV_score <- cnv_score_df[, "CNV_score"]

# ---- CNV density UMAP ----
density_col <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(30))

cnv_density <- plot_density(B_cell_seurat, "CNV_score", size = 0.5,
                            method = "wkde", adjust = 2) +
  scale_color_gradientn(colors = density_col) +
  theme(
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text   = element_blank(),
    aspect.ratio = 1
  )

ggsave(file.path(fig_dir, "fig2d_cnv_density.pdf"),
       egg::set_panel_size(cnv_density, width = unit(16, "cm"), height = unit(16, "cm")),
       width = 24, height = 24, units = "cm")

# ---- CNV ridge (by fine_cell_type) ----
medians <- B_cell_seurat@meta.data %>%
  group_by(fine_cell_type) %>%
  summarise(med = median(CNV_score, na.rm = TRUE), .groups = "drop")

cnv_ridge <- ggplot(B_cell_seurat@meta.data,
                    aes(x = CNV_score, y = 1, fill = fine_cell_type)) +
  geom_density_ridges(alpha = 0.5, scale = 1, size = 0.116) +
  geom_vline(data = medians, aes(xintercept = med, colour = fine_cell_type),
             linetype = "dashed", size = 0.116) +
  scale_fill_manual(values = fine_cell_type_col) +
  scale_color_manual(values = fine_cell_type_col) +
  labs(x = "CNV score", y = NULL) +
  theme_bw() +
  theme(
    text = element_text(size = 5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line  = element_line(colour = "black", size = 0.116),
    axis.ticks = element_blank(),
    axis.text  = element_blank(),
    legend.title = element_text(size = 5)
  ) +
  stat_compare_means(
    comparisons    = list(c("B_Mal", "B_Nor")),
    size           = 5 * 25.4 / 72,
    tip.length     = 0,
    bracket.size   = 0.116,
    method         = "wilcox.test",
    label          = "p.signif"
  )

ggsave(file.path(fig_dir, "fig2d_cnv_ridge.pdf"),
       egg::set_panel_size(cnv_ridge, width = unit(1.5, "cm"), height = unit(1, "cm")),
       width = 10, height = 10, units = "cm")

wilcox_res <- wilcox.test(CNV_score ~ fine_cell_type,
                          data = B_cell_seurat@meta.data, exact = FALSE)
wilcox_res

# ---- AUCell: MZB signatures ----
# reuse B_cell_seurat
count_mat <- as.matrix(B_cell_seurat@assays$RNA@counts)
cells_rankings <- AUCell_buildRankings(count_mat, nCores = 12, plotStats = TRUE)

gene_sets <- list(
  GO_MZB = c("CDH17","DOCK11","PTK2B","DLL1","LFNG","MFNG","NOTCH2","DOCK10","BCL3"),
  Mabbott_MZB = c(
    "ABCB1","RUNDC3B","CYP39A1","PRF1","ADAM28","CDON","MYOF","NEBL",
    "CD1D","ACKR3","FFAR2","S1PR3","ASB2","CREBL2","DPH5","DRD5",
    "DUSP16","GPR156","PIK3R4","PTPN14","TLR3","TSPAN15","ZC3H12C",
    "TERB1","GML","MFHAS1","MS4A7"
  )
)

cells_auc <- AUCell_calcAUC(gene_sets, cells_rankings)
auc_mat   <- as.matrix(getAUC(cells_auc))  # rows = gene sets, cols = cells

# write AUC back to metadata (align by cell names)
auc_df <- t(auc_mat)  # rows = cells, cols = gene sets
auc_df <- auc_df[rownames(B_cell_seurat@meta.data), , drop = FALSE]
B_cell_seurat@meta.data <- cbind(B_cell_seurat@meta.data, auc_df)

# ---- MZB density UMAPs ----
go_mzb_density <- plot_density(B_cell_seurat, "GO_MZB", size = 0.5,
                               method = "wkde", adjust = 2) +
  scale_color_gradientn(colors = density_col) +
  theme(
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text   = element_blank(),
    aspect.ratio = 1
  )

ggsave(file.path(fig_dir, "fig2d_mzb_gobp_density.pdf"),
       egg::set_panel_size(go_mzb_density, width = unit(16, "cm"), height = unit(16, "cm")),
       width = 24, height = 24, units = "cm")

mabbott_mzb_density <- plot_density(B_cell_seurat, "Mabbott_MZB", size = 0.5,
                                    method = "wkde", adjust = 2) +
  scale_color_gradientn(colors = density_col) +
  theme(
    axis.line.x  = element_blank(),
    axis.line.y  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text   = element_blank(),
    aspect.ratio = 1
  )

ggsave(file.path(fig_dir, "fig2d_mzb_mabbott_density.pdf"),
       egg::set_panel_size(mabbott_mzb_density, width = unit(16, "cm"), height = unit(16, "cm")),
       width = 24, height = 24, units = "cm")

# ---- MZB boxplots ----
go_mzb_box <- ggplot(B_cell_seurat@meta.data,
                     aes(x = fine_cell_type, y = GO_MZB, fill = fine_cell_type)) +
  geom_boxplot(outliers = FALSE, size = 0.1, width = 0.5) +
  xlab("") + ylab("GO_MZB score") +
  scale_fill_manual(values = fine_cell_type_col) +
  theme(
    text = element_text(size = 5),
    panel.grid.major = element_line(colour = "grey90", size = 0.116),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line  = element_line(colour = "black", size = 0.116),
    axis.ticks = element_line(colour = "black", size = 0.116),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(colour = "black", angle = 60, vjust = 1, hjust = 1, size = 5),
    axis.text.y  = element_blank(),
    legend.position = "none",
    legend.title    = element_text(size = 5)
  ) +
  stat_compare_means(
    comparisons  = list(c("B_Mal", "B_Nor")),
    size         = 5 * 25.4 / 72,
    tip.length   = 0,
    bracket.size = 0.116,
    method       = "wilcox.test",
    label        = "p.signif"
  ) +
  coord_cartesian(clip = "off")

ggsave(file.path(fig_dir, "fig2d_mzb_gobp_boxplot.pdf"),
       egg::set_panel_size(go_mzb_box, width = unit(1, "cm"), height = unit(1, "cm")),
       width = 6, height = 6, units = "cm")

mabbott_mzb_box <- ggplot(B_cell_seurat@meta.data,
                          aes(x = fine_cell_type, y = Mabbott_MZB, fill = fine_cell_type)) +
  geom_boxplot(outliers = FALSE, size = 0.1, width = 0.4) +
  xlab("") + ylab("Mabbott_MZB score") +
  scale_fill_manual(values = fine_cell_type_col) +
  theme(
    text = element_text(size = 5),
    panel.grid.major = element_line(colour = "grey90", size = 0.116),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line  = element_line(colour = "black", size = 0.116),
    axis.ticks = element_line(colour = "black", size = 0.116),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(colour = "black", angle = 60, vjust = 1, hjust = 1, size = 5),
    axis.text.y  = element_blank(),
    legend.position = "none",
    legend.title   = element_text(size = 5)
  ) +
  stat_compare_means(
    comparisons  = list(c("B_Mal", "B_Nor")),
    size         = 5 * 25.4 / 72,
    tip.length   = 0,
    bracket.size = 0.116,
    method       = "wilcox.test",
    label        = "p.signif"
  ) +
  coord_cartesian(clip = "off")

ggsave(file.path(fig_dir, "fig2d_mzb_mabbott_boxplot.pdf"),
       egg::set_panel_size(mabbott_mzb_box, width = unit(1, "cm"), height = unit(1, "cm")),
       width = 6, height = 6, units = "cm")