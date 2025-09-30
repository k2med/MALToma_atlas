# ---- packages ----
library(edgeR)
library(dplyr)
library(clusterProfiler)
library(DOSE)
library(msigdbr)
library(ggplot2)
library(egg)

# ---- paths & sources ----
src_dir <- "../source"
fig_dir <- "."

# ---- inputs ----
genesets <- readRDS(file.path(src_dir, "msigdbr_gene_sets.rds"))

genelist_df <- read.csv("MP_list.csv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
genelist <- as.list(genelist_df)  # each column becomes a vector

use_genesets <- c(
  "Dna replication [GOBP]",
  "Mitotic g1 phase and g1 s transition [REACTOME]",
  "Endoplasmic reticulum protein containing complex [GOCC]",
  "Unfolded protein response upr [REACTOME]",
  "Chromosome segregation [GOBP]",
  "Mitotic spindle [HALLMARK]",
  "B cell receptor signaling pathway [GOBP]",
  "B cell mediated immunity [GOBP]",
  "Spliceosomal complex [GOCC]",
  "Ribonucleoprotein complex binding [GOMF]",
  "Immunoglobulin complex circulating [GOCC]",
  "Immunoglobulin receptor binding [GOMF]",
  "Il2 stat5 signaling [HALLMARK]",
  "Chemokine receptors bind chemokines [REACTOME]"
)

# subset pathways and build TERM2GENE with required two columns
genesets_sub <- subset(genesets, pathID_rename %in% use_genesets)
genesets_sub$pathID <- genesets_sub$pathID_rename

# ---- enrichment ----
enrich_genesets <- function(term2gene, gene_vec) {
  er <- enricher(
    gene          = gene_vec,
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    minGSSize     = 10,
    maxGSSize     = 500,
    qvalueCutoff  = 0.20,
    TERM2GENE     = term2gene,
    TERM2NAME     = NA
  )
  if (is.null(er)) return(NULL)
  er@result
}

enrich_result <- lapply(genelist, function(gv) {
  enrich_genesets(genesets_sub, gv)
})

# ---- combine results ----
results_df <- bind_rows(
  lapply(names(enrich_result), function(nm) {
    df <- enrich_result[[nm]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$Group <- nm
    df
  })
)

results_df$ID <- factor(results_df$ID, levels = rev(use_genesets))

# ---- plot: barplot ----
enrich_barplot <- ggplot(data = results_df, aes(x = -log10(pvalue), y = ID)) +
  geom_col(fill = "#7f9db4", width = 0.5) +
  facet_grid(. ~ Group, scales = "free") +
  xlab("-log10(P value)") + ylab(NULL) +
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
    axis.text.y = element_text(colour = "black", size = 5),
    legend.position = "right",
    legend.title = element_text(size = 5),
    strip.background = element_blank(),
    strip.text = element_text(size = 5)
  )

ggsave(
  filename = file.path(fig_dir, "fig3b.pdf"),
  plot = egg::set_panel_size(
    enrich_barplot,
    width  = grid::unit(0.6, "cm"),
    height = grid::unit(3.0, "cm")
  ),
  width = 16, height = 10, units = "cm"
)
