# ---- packages ----
library(Seurat)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(egg)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- homologene mapping (human-mouse orthologs) ----
homologene_file <- file.path(src_dir, "homologene.data")
columns <- c("HomoloGeneID","TaxonomyID","GeneID","GeneSymbol","GeneDescription","ProteinGI","ProteinAccession")
homologene_df <- read_tsv(homologene_file, col_names = columns, col_types = cols())

mapping_df <- homologene_df %>%
  filter(TaxonomyID %in% c(9606, 10090)) %>%
  group_by(HomoloGeneID) %>%
  filter(all(c(9606, 10090) %in% TaxonomyID)) %>%
  summarise(
    HumanGeneSymbol = first(GeneSymbol[TaxonomyID == 9606]),
    MouseGeneSymbol = first(GeneSymbol[TaxonomyID == 10090]),
    HumanGeneID     = first(GeneID[TaxonomyID == 9606]),
    MouseGeneID     = first(GeneID[TaxonomyID == 10090]),
    .groups = "drop"
  )

# ---- inputs ----
seurat_myeloid <- readRDS(file.path(data_dir, "snrna_myeloid_subset_seurat.rds"))

# ---- helpers ----
# 1) DEG per subpopulation set, map to mouse, export lists for IREA website
deg_and_export_for_irea <- function(seurat_obj, subtypes_vec, out_prefix,
                                    min_pct = 0.25, logfc_thr = 0.25) {
  obj <- subset(seurat_obj, fine_cell_type %in% subtypes_vec)
  obj$fine_cell_type <- droplevels(obj$fine_cell_type)
  Idents(obj) <- obj$fine_cell_type
  
  deg_df <- FindAllMarkers(
    obj, assay = "RNA", only.pos = FALSE, test.use = "MAST",
    min.pct = min_pct, logfc.threshold = logfc_thr
  )
  
  # merge human -> mouse symbol
  deg_df <- deg_df %>%
    left_join(mapping_df, by = c("gene" = "HumanGeneSymbol")) %>%
    mutate(mouse_gene = if_else(is.na(MouseGeneSymbol), "", MouseGeneSymbol))
  
  # filter DEGs
  deg_filtered <- deg_df %>%
    filter(avg_log2FC > logfc_thr,
           (pct.1 > min_pct | pct.2 > min_pct),
           p_val_adj < 0.05)
  
  # ---- IMPORTANT: website step (IREA) ----
  # Export gene lists for each subtype to upload to https://www.immune-dictionary.org
  # Perform the hypergeometric test on the website.
  # After running on the website, download CSV results named "<LABEL>_IREA_output.csv"
  # where <LABEL> is a concise subtype label.
  for (ct in levels(obj$fine_cell_type)) {
    mouse_ct <- deg_filtered %>% filter(cluster == ct, mouse_gene != "") %>% pull(mouse_gene) %>% unique()
    write_lines(mouse_ct, file.path(fig_dir, paste0(out_prefix, "_", ct, "_DEG.txt")))
  }
  
  return(invisible(NULL))
}

# 2) read IREA CSVs, clean greek letters, filter cytokines, scale ES, compute -log10(padj)
replace_greek_to_ascii <- function(x) {
  greek <- "αβγ"
  ascii <- "abg"
  chartr(greek, ascii, x)
}

read_irea_and_prepare <- function(file_map, use_cytokines, cell_type_levels) {
  # file_map: named character vector, names = display cell type, values = CSV filenames
  bind_rows(lapply(names(file_map), function(ct) {
    df <- read.csv(file_map[[ct]], fileEncoding = "UTF-8", stringsAsFactors = FALSE)
    df$Cytokine <- replace_greek_to_ascii(df$Cytokine)
    df <- df[df$Cytokine %in% use_cytokines, , drop = FALSE]
    df$ES_scaled <- as.numeric(scale(df$Enrichment.Score))
    df$log10padj <- -log10(df$padj)
    df$cell_type <- ct
    df
  })) %>%
    mutate(
      Cytokine = factor(Cytokine, levels = rev(use_cytokines)),
      cell_type = factor(cell_type, levels = cell_type_levels)
    )
}

# 3) dotplot for IREA results
plot_irea_dot <- function(df, size_limits = c(0, 30), outfile, fig_dir) {
  p <- ggplot(df, aes(x = cell_type, y = Cytokine)) +
    geom_point(aes(color = ES_scaled, size = log10padj), shape = 19, stroke = 0.01) +
    labs(x = NULL, y = NULL) +
    scale_color_gradientn(colors = rev(hcl.colors(n = 10, palette = "Earth"))) +
    scale_size(range = c(0, 2.4), limits = size_limits,
               breaks = pretty(size_limits, n = 4)) +
    theme_bw(base_size = 5) +
    theme(
      text = element_text(size = 5),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.116),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", linewidth = 0.116),
      axis.line = element_line(colour = "black", linewidth = 0.116),
      axis.ticks = element_line(colour = "black", linewidth = 0.116),
      axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust = 1, size = 5),
      axis.text.y = element_text(colour = "black", size = 5),
      legend.position = "right",
      legend.title = element_text(size = 5)
    ) +
    guides(
      size = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(stroke = 0.4)),
      fill = guide_colourbar(title.position = "top", title.hjust = 0.5)
    )
  
  ggsave(
    filename = file.path(fig_dir, outfile),
    plot = egg::set_panel_size(
      p,
      width  = grid::unit(length(unique(df$cell_type)) / 5, "cm"),
      height = grid::unit(length(levels(df$Cytokine)) / 5, "cm")
    ),
    width = 10, height = 10, units = "cm"
  )
}

# ========== DC =============
dc_subtypes <- c("cDC_XCR1", "cDC_CLEC10A", "pDC_CLEC4C")
deg_and_export_for_irea(seurat_myeloid, dc_subtypes, out_prefix = "DC")

# ---- WEBSITE STEP (do this outside R) -------------------------
# Open https://www.immune-dictionary.org
# For each DC subtype (cDC_XCR1, cDC_CLEC10A, pDC_CLEC4C),
# upload the corresponding DEG list file exported above:
#   DC_cDC_XCR1_DEG.txt
#   DC_cDC_CLEC10A_DEG.txt
#   DC_pDC_CLEC4C_DEG.txt
# Run the hypergeometric test on the website.
# Download CSV results named, for example:
#   XCR1_IREA_output.csv
#   CLEC10A_IREA_output.csv
#   CLEC4C_IREA_output.csv
# ---------------------------------------------------------------

# After website processing, read CSVs and plot
dc_files <- c(
  "cDC_XCR1"   = "XCR1_IREA_output.csv",
  "cDC_CLEC10A"= "CLEC10A_IREA_output.csv",
  "pDC_CLEC4C" = "CLEC4C_IREA_output.csv"
)

dc_use_cytokines <- c("IFN-a1","IFN-b","IFN-g","IL-1a","IL-1b","TNF-a","GM-CSF")
dc_df <- read_irea_and_prepare(file_map = dc_files,
                               use_cytokines = dc_use_cytokines,
                               cell_type_levels = names(dc_files))
plot_irea_dot(dc_df, size_limits = c(0, 30), outfile = "fig4k.pdf", fig_dir = fig_dir)

# ======== Macrophage ==========
mac_subtypes <- c("Mac_SPARC", "Mac_STAB1", "Mac_MARCO", "Mac_CHIT1")
deg_and_export_for_irea(seurat_myeloid, mac_subtypes, out_prefix = "Mac")

# ---- WEBSITE STEP (do this outside R) -------------------------
# Open https://www.immune-dictionary.org
# For each Macrophage subtype (Mac_SPARC, Mac_STAB1, Mac_MARCO, Mac_CHIT1),
# upload the corresponding DEG file exported above:
#   Mac_Mac_SPARC_DEG.txt
#   Mac_Mac_STAB1_DEG.txt
#   Mac_Mac_MARCO_DEG.txt
#   Mac_Mac_CHIT1_DEG.txt
# Run the hypergeometric test on the website.
# Download CSV results named, for example:
#   SPARC_IREA_output.csv
#   STAB1_IREA_output.csv
#   MARCO_IREA_output.csv
#   CHIT1_IREA_output.csv
# ---------------------------------------------------------------

mac_files <- c(
  "Mac_SPARC" = "SPARC_IREA_output.csv",
  "Mac_STAB1" = "STAB1_IREA_output.csv",
  "Mac_MARCO" = "MARCO_IREA_output.csv",
  "Mac_CHIT1" = "CHIT1_IREA_output.csv"
)

mac_use_cytokines <- c("IL-15","TNF-a","IL-1a","IL-1b","IL-4","IL-13","M-CSF")
mac_df <- read_irea_and_prepare(file_map = mac_files,
                                use_cytokines = mac_use_cytokines,
                                cell_type_levels = names(mac_files))
plot_irea_dot(mac_df, size_limits = c(0, 15), outfile = "figS4b.pdf", fig_dir = fig_dir)