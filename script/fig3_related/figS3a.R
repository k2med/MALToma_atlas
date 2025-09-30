# ---- packages ----
library(Seurat)
library(msigdbr)
library(tidyverse)
library(clusterProfiler)
library(dplyr)
library(enrichplot)
library(GseaVis)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir <- "../source"
fig_dir <- "."

source(file.path(src_dir, "custom_colors.R"))
source(file.path(src_dir, "plot_gsea.R"))

# ---- inputs ----
B_cell_seurat <- readRDS(file.path(data_dir, "snrna_b_subset_seurat.rds"))

B_cell_seurat@active.ident <- factor(B_cell_seurat$fine_cell_type)

fine_cell_type_markers <- FindMarkers(B_cell_seurat, ident.1 = 'B_Mal', ident.2 = 'B_Nor',
                                      assay = 'RNA',
                                      only.pos = F, 
                                      test.use = 'MAST',
                                      min.pct = 0,
                                      logfc.threshold = 0, 
                                      min.diff.pct = 0)

gene_logFC <- fine_cell_type_markers$avg_log2FC
names(gene_logFC) <- rownames(fine_cell_type_markers)
sort_gene_logFC <- gene_logFC[order(gene_logFC, decreasing = T)]

genesets <- msigdbr(species = "Homo sapiens"
                    , category = "H") %>%
  subset(select = c("gs_name","gene_symbol")) %>%
  as.data.frame() %>% unique()

gsea <- GSEA(sort_gene_logFC,
             TERM2GENE = genesets,
             pvalueCutoff = 1,
             pAdjustMethod = 'BH')

for (genelist.i in c("HALLMARK_G2M_CHECKPOINT","HALLMARK_E2F_TARGETS","HALLMARK_MITOTIC_SPINDLE")) {
  
  gseaNb(
    object= gsea,
    geneSetID= genelist.i,
    subPlot= 2,
    addPval= T,
    htHeight=0.3,
    pvalX= 0.95,
    pvalY= 0.8
  )
  ggsave(paste0('figS3a_', genelist.i, '.pdf'), height = 4, width = 6)
}


genesets <- msigdbr(species = "Homo sapiens"
                    , category = "C2", subcategory = "KEGG") %>%
  subset(select = c("gs_name","gene_symbol")) %>%
  as.data.frame() %>% unique()

gsea <- GSEA(sort_gene_logFC,
             TERM2GENE = genesets,
             pvalueCutoff = 1,
             pAdjustMethod = 'BH')

for (genelist.i in c("KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY")) {
  
  gseaNb(
    object= gsea,
    geneSetID= genelist.i,
    subPlot= 2,
    addPval= T,
    htHeight=0.3,
    pvalX= 0.95,
    pvalY= 0.8
  )
  ggsave(paste0('figS3a_', genelist.i, '.pdf'), height = 4, width = 6)
}
