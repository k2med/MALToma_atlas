# ---- packages ----
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(egg)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))
source(file.path(src_dir, "plot_heatmap.R"))

# ---- function ----
run_marker_heatmap <- function(
    seurat_rds_path,
    celltype_levels,
    column_title,
    outfile_pdf,
    assay_name        = "RNA",
    test_method       = "MAST",
    min_pct           = 0.1,
    logfc_threshold   = 0.1,
    padj_thresh       = 0.05,
    log2fc_sig        = 0.25,
    top_n_left_labels = 5,
    hm_name           = "Markers_heatmap",
    hm_legend_name    = "Relative expression",
    hm_palette        = heatmap_col3
) {
  # load object
  obj <- readRDS(seurat_rds_path)
  
  # enforce factor order for fine_cell_type
  obj$fine_cell_type <- factor(obj$fine_cell_type, levels = celltype_levels)
  Idents(obj) <- obj$fine_cell_type
  
  # find markers
  markers_all <- FindAllMarkers(
    object = obj,
    assay = assay_name,
    only.pos = FALSE,
    test.use = test_method,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold
  )
  
  # significant markers and unique best cluster per gene
  markers_sig <- subset(markers_all, (p_val_adj < padj_thresh) & (avg_log2FC > log2fc_sig))
  
  markers_sig_unique <- markers_sig %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(avg_log2FC, n = 1) %>%
    dplyr::ungroup()
  
  ord_df <- as.data.frame(
    markers_sig_unique[order(markers_sig_unique$cluster,
                             markers_sig_unique$avg_log2FC,
                             decreasing = c(FALSE, TRUE)), ]
  )
  rownames(ord_df) <- ord_df$gene
  
  # shuffle columns within each group to avoid long runs of the same type
  meta <- obj@meta.data
  meta_shuffled <- meta %>%
    dplyr::mutate(row_id = rownames(meta)) %>%
    dplyr::group_by(fine_cell_type) %>%
    dplyr::group_modify(~ dplyr::slice_sample(.x, n = nrow(.x))) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  rownames(meta_shuffled) <- meta_shuffled$row_id
  
  # expression matrix aligned: genes × cells
  expr_mat <- as.matrix(
    obj@assays[[assay_name]]@data[rownames(ord_df), rownames(meta_shuffled), drop = FALSE]
  )
  
  # annotations
  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Cell_type = meta_shuffled$fine_cell_type,
    col = list(Cell_type = fine_cell_type_col),
    annotation_name_gp = grid::gpar(fontface = "plain")
  )
  
  # heatmap object (right)
  plot_markers_heatmap <- heatmap_simple(
    expr_mat,
    top_annotation = top_annotation,
    column_title   = column_title,
    name           = hm_name,
    scale_rows     = TRUE,
    width          = grid::unit(3, "in"),
    height         = grid::unit(6, "in"),
    raster_quality = 7,
    legend_name    = hm_legend_name,
    color_palette  = hm_palette,
    show_annotation_name = FALSE,
    color_range    = 2 * seq(-1, 1, length.out = length(hm_palette) - 1)
  )
  
  # draw: left labels + heatmap
  pdf(outfile_pdf, width = 10, height = 12)
  suppressWarnings({
    ComplexHeatmap::draw(
      top_markers_left(expr_mat, top_n = top_n_left_labels,
                       ord = ord_df, ord_col = "cluster") + plot_markers_heatmap,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      merge_legends = FALSE
    )
  })
  
  # decorate rectangles
  rect_coords <- rectangle_annotation_coordinates(
    ord_df$cluster,
    meta_shuffled$fine_cell_type
  )
  
  ComplexHeatmap::decorate_heatmap_body(hm_name, {
    grid::grid.rect(
      x = grid::unit(rect_coords$x, "native"),
      y = grid::unit(rect_coords$y, "native"),
      width  = grid::unit(rect_coords$w, "native"),
      height = grid::unit(rect_coords$h, "native"),
      hjust = 0, vjust = 1,
      gp = grid::gpar(col = "white", lty = 1, lwd = 1)
    )
  })
  dev.off()
}

# ---- calls ----
run_marker_heatmap(
  seurat_rds_path = file.path(data_dir, "snrna_t_subset_seurat.rds"),
  celltype_levels = c("CD8T_CCL5","CD8T_GNLY","CD8T_HAVCR2","CD4T_FOXP3","CD4T_CXCL13","CD4T_IL7R"),
  column_title    = "T cells",
  outfile_pdf     = file.path(fig_dir, "figS2c_t_cell.pdf")
)

run_marker_heatmap(
  seurat_rds_path = file.path(data_dir, "snrna_myeloid_subset_seurat.rds"),
  celltype_levels = c("Mac_SPARC","Mac_STAB1","Mac_MARCO","Mac_CHIT1",
                      "cDC_XCR1","cDC_CLEC10A","pDC_CLEC4C","Neu_G0S2","Mast_CPA3"),
  column_title    = "Myeloid cells",
  outfile_pdf     = file.path(fig_dir, "figS2c_myeloid_cell.pdf")
)

run_marker_heatmap(
  seurat_rds_path = file.path(data_dir, "snrna_stromal_subset_seurat.rds"),
  celltype_levels = c("Endo_ACKR1","Endo_CXCL12","Endo_CCL21","Endo_HPGD","Peri_PDGFRB",
                      "Fibro_PI16","Fibro_CCL19","Fibro_LRRC15","Fibro_TCF21","FDC_CR2"),
  column_title    = "Stromal cells",
  outfile_pdf     = file.path(fig_dir, "figS2c_stromal_cell.pdf")
)

run_marker_heatmap(
  seurat_rds_path = file.path(data_dir, "snrna_alveolar_subset_seurat.rds"),
  celltype_levels = c("AT1_AGER","AT2_SFTPC"),
  column_title    = "Alveolar cells",
  outfile_pdf     = file.path(fig_dir, "figS2c_alveolar_cell.pdf")
)