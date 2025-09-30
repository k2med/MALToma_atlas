# ---- packages ----
library(Seurat)
library(SeuratObject)
library(nichenetr)
library(tidyverse)
library(ggplot2)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

# ---- load nichenet resources (human/mouse) ----
# Files for v2: https://zenodo.org/records/7074291
organism <- "human"
if (organism == "human") {
  lr_network <- readRDS("./lr_network_human_21122021.rds")
  ligand_target_matrix <- readRDS("./ligand_target_matrix_nsga2r_final.rds")
  weighted_networks <- readRDS("./weighted_networks_nsga2r_final.rds")
} else if (organism == "mouse") {
  lr_network <- readRDS("./lr_network_mouse_21122021.rds")
  ligand_target_matrix <- readRDS("./ligand_target_matrix_nsga2r_final_mouse.rds")
  weighted_networks <- readRDS("./weighted_networks_nsga2r_final_mouse.rds")
} else {
  stop("Unsupported organism; choose 'human' or 'mouse'.")
}
lr_network <- lr_network %>% distinct(from, to)

# ---- load Seurat all-cells object & subset usable cells ----
seurat_all_cells <- readRDS(file.path(data_dir, "snrna_all_cells_seurat.rds"))
seurat_use_cells <- subset(
  seurat_all_cells,
  fine_cell_type %in% c("Fibro_TCF21", "AT1_AGER", "AT2_SFTPC", "Hepa_ALB", "Endo_HPGD"),
  invert = TRUE
)
seurat_use_cells$fine_cell_type <- droplevels(seurat_use_cells$fine_cell_type)
Idents(seurat_use_cells) <- seurat_use_cells$fine_cell_type

# ---- helper: expressed genes & potential ligands ----
get_receiver_context <- function(receiver_types, seurat_obj = seurat_use_cells, pct = 0.10) {
  expressed_genes <- nichenetr::get_expressed_genes(
    receiver_types,   # celltype_oi
    seurat_obj,       # Seurat object
    pct = pct
  )
  
  all_receptors <- unique(lr_network$to)
  expressed_receptors <- intersect(all_receptors, expressed_genes)
  potential_ligands <- lr_network %>%
    dplyr::filter(to %in% expressed_receptors) %>%
    dplyr::pull(from) %>% unique()
  
  list(
    expressed_genes     = expressed_genes,
    potential_ligands   = potential_ligands,
    seurat_obj_receiver = subset(seurat_obj, idents = receiver_types)
  )
}

# ---- helper: differential geneset in receiver (ident_1 vs ident_2) ----
get_receiver_geneset <- function(seurat_obj_receiver, ident_1, ident_2,
                                 fdr_cutoff = 0.05, lfc_cutoff = 0.25) {
  de_tbl <- FindMarkers(
    object = seurat_obj_receiver,
    ident.1 = ident_1,
    ident.2 = ident_2,
    group.by = "fine_cell_type",
    test.use = "MAST",
    logfc.threshold = 0,
    min.pct = 0
  ) %>% rownames_to_column("gene")
  
  geneset_oi <- de_tbl %>%
    filter(p_val_adj <= fdr_cutoff, avg_log2FC >= lfc_cutoff) %>%
    pull(gene) %>%
    intersect(rownames(ligand_target_matrix))
  
  geneset_oi
}

# ---- helper: ligand activities (AUPR) ----
rank_ligands <- function(geneset_oi, expressed_genes_receiver, potential_ligands) {
  if (length(geneset_oi) == 0) return(NULL)
  
  background_expressed <- intersect(expressed_genes_receiver, rownames(ligand_target_matrix))
  
  ligand_activities <- predict_ligand_activities(
    geneset                      = geneset_oi,
    background_expressed_genes   = background_expressed,
    ligand_target_matrix         = ligand_target_matrix,
    potential_ligands            = potential_ligands
  ) %>%
    arrange(desc(aupr_corrected)) %>%
    mutate(rank = rank(-aupr_corrected, ties.method = "min"))
  
  ligand_activities
}

# ---- helper: plot ligand AUPR heatmap (left panel) ----
plot_ligand_aupr <- function(ligand_activities, top_k, out_pdf) {
  if (is.null(ligand_activities) || nrow(ligand_activities) == 0) {
    warning("No ligand_activities available; skip AUPR plot.")
    return(invisible(NULL))
  }
  best_ligands <- ligand_activities %>%
    slice_max(order_by = aupr_corrected, n = top_k, with_ties = FALSE) %>%
    arrange(desc(aupr_corrected)) %>%
    pull(test_ligand)
  
  vis_ligand_aupr <- ligand_activities %>%
    filter(test_ligand %in% best_ligands) %>%
    column_to_rownames("test_ligand") %>%
    select(aupr_corrected) %>%
    arrange(aupr_corrected) %>%
    as.matrix(ncol = 1)
  
  p <- make_heatmap_ggplot(
    vis_ligand_aupr,
    "Prioritized ligands", "Ligand activity",
    legend_title = "AUPR", color = "darkorange"
  ) + theme(axis.text.x.top = element_blank())
  
  ggsave(file.path(fig_dir, out_pdf), p, width = 6, height = 12, units = "cm")
  invisible(best_ligands)
}

# ---- helper: plot ligand→target heatmap (right panel) ----
plot_ligand_targets <- function(best_ligands, geneset_oi, out_pdf,
                                n_targets_per_ligand = 250,
                                rp_cutoff = 0.25,
                                target_n_first = 20,       # set to NULL to keep all
                                low_col = "whitesmoke",
                                high_col = "purple") {
  if (length(best_ligands) == 0 || length(geneset_oi) == 0) {
    warning("Empty best_ligands or geneset_oi; skip target plot.")
    return(invisible(NULL))
  }
  
  links_df <- best_ligands %>%
    lapply(
      get_weighted_ligand_target_links,
      geneset               = geneset_oi,
      ligand_target_matrix  = ligand_target_matrix,
      n                     = n_targets_per_ligand
    ) %>%
    bind_rows() %>%
    drop_na()
  
  links_mat <- prepare_ligand_target_visualization(
    ligand_target_df        = links_df,
    ligand_target_matrix    = ligand_target_matrix,
    cutoff                  = rp_cutoff
  )
  
  order_lig <- intersect(best_ligands, colnames(links_mat)) %>% rev()
  order_tar <- links_df$target %>% unique() %>% intersect(rownames(links_mat))
  
  if (length(order_lig) == 0 || length(order_tar) == 0) {
    warning("No ligand-target entries after filtering; skip target plot.")
    return(invisible(NULL))
  }
  
  if (!is.null(target_n_first) && target_n_first > 0 && length(order_tar) > target_n_first) {
    order_tar <- order_tar[1:target_n_first]
  }
  
  vis_mat <- t(links_mat[order_tar, order_lig, drop = FALSE])
  
  p <- make_heatmap_ggplot(
    vis_mat,
    "Prioritized ligands", "Predicted target genes",
    color = high_col, legend_title = "Regulatory potential"
  ) + scale_fill_gradient2(low = low_col, high = high_col)
  
  ggsave(file.path(fig_dir, out_pdf), p, width = 24, height = 12, units = "cm")
  invisible(TRUE)
}

# ---- wrapper: one receiver contrast (calls the short helpers) ----
run_receiver_contrast <- function(receiver_types,
                                  ident_1, ident_2,
                                  top_k,
                                  out_prefix,
                                  target_n_first = 20,
                                  rp_cutoff = 0.25,
                                  seurat_obj = seurat_use_cells) {
  ctx <- get_receiver_context(receiver_types, seurat_obj = seurat_obj, pct = 0.10)
  
  geneset_oi <- get_receiver_geneset(
    seurat_obj_receiver = ctx$seurat_obj_receiver,
    ident_1 = ident_1,
    ident_2 = ident_2,
    fdr_cutoff = 0.05,
    lfc_cutoff = 0.25
  )
  if (length(geneset_oi) == 0) {
    warning(sprintf("No geneset_oi for %s vs %s; skip.", ident_1, paste(ident_2, collapse = ",")))
    return(invisible(NULL))
  }
  
  ligand_activities <- rank_ligands(
    geneset_oi               = geneset_oi,
    expressed_genes_receiver = ctx$expressed_genes,
    potential_ligands        = ctx$potential_ligands
  )
  
  best_ligands <- plot_ligand_aupr(
    ligand_activities = ligand_activities,
    top_k             = top_k,
    out_pdf           = paste0(out_prefix, "_left.pdf")
  )
  if (is.null(best_ligands) || length(best_ligands) == 0) return(invisible(NULL))
  
  plot_ligand_targets(
    best_ligands        = best_ligands,
    geneset_oi          = geneset_oi,
    out_pdf             = paste0(out_prefix, "_right.pdf"),
    n_targets_per_ligand = 250,
    rp_cutoff           = rp_cutoff,
    target_n_first      = target_n_first
  )
  invisible(TRUE)
}

# ---- block A: Fibro_CCL19 vs (Fibro_PI16, Fibro_LRRC15) -> fig5c_* ----
run_receiver_contrast(
  seurat_obj      = seurat_use_cells,
  receiver_types  = c("Fibro_PI16", "Fibro_CCL19", "Fibro_LRRC15"),
  ident_1         = "Fibro_CCL19",
  ident_2         = c("Fibro_PI16", "Fibro_LRRC15"),
  top_k           = 5,
  out_prefix      = "fig5c",
  target_n_first  = 20,
  rp_cutoff       = 0.25
)

# ---- block B: Fibro_LRRC15 vs (Fibro_PI16, Fibro_CCL19) -> fig5d_* ----
run_receiver_contrast(
  receiver_types  = c("Fibro_PI16", "Fibro_CCL19", "Fibro_LRRC15"),
  ident_1         = "Fibro_LRRC15",
  ident_2         = c("Fibro_PI16", "Fibro_CCL19"),
  top_k           = 5,
  out_prefix      = "fig5d",
  target_n_first  = 20,
  rp_cutoff       = 0.25
)

# ---- block C: Endo_CXCL12 vs Endo_ACKR1 -> fig5i_* ----
run_receiver_contrast(
  receiver_types  = c("Endo_ACKR1", "Endo_CXCL12"),
  ident_1         = "Endo_CXCL12",
  ident_2         = "Endo_ACKR1",
  top_k           = 7,
  out_prefix      = "fig5i",
  target_n_first  = NULL,    # keep all predicted targets in the right panel
  rp_cutoff       = 0.25
)