# ---- packages ----
library(tidyverse)
library(Seurat)
library(mistyR)
library(progeny)
library(argparse)

# ---- paths & sources ----
data_dir <- "../../data/spatial_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))
source(file.path(src_dir, "run_misty_seurat.R"))

# ---- discover samples ----
sample_ids <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
slide_ids   <- sample_ids

# ---- helper: read Visium image ----
read_visium_image_custom <- function(image_dir,
                                     image_name   = "tissue_lowres_image.png",
                                     filter_in_tissue = TRUE) {
  img            <- readPNG(source = file.path(image_dir, image_name))
  scale_factors  <- fromJSON(txt = file.path(image_dir, "scalefactors_json.json"))
  pos_path       <- Sys.glob(paths = file.path(image_dir, "tissue_positions*"))
  pos_df <- read.csv(
    file      = pos_path,
    col.names = c("barcode", "in_tissue", "row", "col", "imagerow", "imagecol"),
    header    = ifelse(basename(pos_path) == "tissue_positions.csv", TRUE, FALSE),
    as.is     = TRUE,
    row.names = 1
  )
  if (filter_in_tissue) {
    pos_df <- pos_df[pos_df$in_tissue == 1, , drop = FALSE]
  }
  # compute spot radius on the plotted image
  unnorm_radius <- scale_factors$fiducial_diameter_fullres * scale_factors$tissue_lowres_scalef
  spot_radius   <- unnorm_radius / max(dim(img))
  
  return(
    new(
      Class         = "VisiumV1",
      image         = img,
      scale.factors = scalefactors(
        spot     = scale_factors$spot_diameter_fullres,
        fiducial = scale_factors$fiducial_diameter_fullres,
        hires    = scale_factors$tissue_hires_scalef,
        lowres   = scale_factors$tissue_lowres_scalef
      ),
      coordinates   = pos_df,
      spot.radius   = spot_radius
    )
  )
}

# ---- load all sections into a single Seurat object ----
seurat_by_sample <- list()

for (sid in sample_ids) {
  sample_dir <- file.path(data_dir, sid)
  
  expr_mat <- Read10X_h5(file.path(sample_dir, "filtered_feature_bc_matrix.h5"))
  
  img <- read_visium_image_custom(
    image_dir       = file.path(sample_dir, "spatial"),
    filter_in_tissue= TRUE,
    image_name      = "tissue_hires_image.png"
  )
  # subset image spots to match barcodes in the matrix
  img <- img[colnames(expr_mat)]
  DefaultAssay(img) <- "Spatial"
  # align hires/lowres scale factors to avoid downstream issues
  img@scale.factors$lowres <- img@scale.factors$hires
  
  # prefix barcodes/coordinates with sample id to keep them unique after merge
  rownames(img@coordinates) <- paste(sid, rownames(img@coordinates), sep = "_")
  colnames(expr_mat)        <- paste(sid, colnames(expr_mat),      sep = "_")
  
  seu_i <- CreateSeuratObject(
    counts      = expr_mat,
    assay       = "Spatial",
    names.field = 1:2,
    names.delim = "_"
  )
  seu_i[[sid]] <- img
  
  seurat_by_sample[[sid]] <- seu_i
}

seurat_spatial <- merge(seurat_by_sample[[1]], y = seurat_by_sample[-1])

# ---- load cell2location results (spot x celltype abundances) ----
cell2loc_df <- read.csv("./cell2location_map/st_cell2location_res.csv",
                        row.names = 1, check.names = FALSE)

# clean column names (remove prefixes/special chars)
colnames(cell2loc_df) <- stringr::str_replace_all(
  colnames(cell2loc_df),
  c("q05cell_abundance_w_sf_" = "", "\\+" = "", " " = "-")
)

# ---- helper for MISTy colocalization ----
run_colocalization <- function(
    slide,
    assay,
    useful_features,
    out_label,
    misty_out_alias = "misty_out"
) {
  # define views
  view_assays   <- list(main = assay, juxta = assay, para = assay)
  view_features <- list(main = useful_features, juxta = useful_features, para = useful_features)
  view_types    <- list(main = "intra", juxta = "juxta", para = "para")
  # extra params: l for para, n for juxta
  view_params   <- list(main = NULL, juxta = 2, para = 5)
  
  out_prefix <- paste0(misty_out_alias, out_label, "_", assay)
  
  run_misty_seurat(
    visium.slide  = slide,
    view.assays   = view_assays,
    view.features = view_features,
    view.types    = view_types,
    view.params   = view_params,
    spot.ids      = NULL,
    out.alias     = out_prefix
  )
  
  return(out_prefix)
}

# ---- figS7d top (ligand-cell abundance) ----
out_dir_top <- "figS7d_top_"

cellchat_db   <- CellChatDB.human$interaction
slide_markers <- unique(cellchat_db[cellchat_db$pathway_name %in% c("CD40", "BAFF"), "ligand"])

for (slide_id in slide_ids) {
  slide <- subset(seurat_spatial, orig.ident == slide_id)
  slide <- NormalizeData(slide, verbose = FALSE)
  
  cell_mtx   <- t(cell2loc_df[colnames(slide), ])
  expr_mat   <- as.matrix(GetAssayData(slide, slot = "data", assay = "Spatial"))
  ligand_mtx <- as.matrix(expr_mat[slide_markers, , drop = FALSE])
  merge_mtx  <- rbind(cell_mtx, ligand_mtx)
  
  slide@assays$predictions <- CreateAssayObject(counts = merge_mtx)
  DefaultAssay(slide)      <- "predictions"
  slide@images             <- seurat_spatial@images[slide_id]
  
  useful_features <- rownames(slide)
  
  invisible(run_colocalization(
    slide            = slide,
    useful_features  = useful_features,
    out_label        = slide_id,
    assay            = "predictions",
    misty_out_alias  = out_dir_top
  ))
}

# collect results
misty_res_lme1 <- collect_results(c(
  paste0("./", out_dir_top, "S31_ST_predictions/"),
  paste0("./", out_dir_top, "S32_ST_predictions/"),
  paste0("./", out_dir_top, "S71_ST_predictions/")
))
misty_res_lme2 <- collect_results(c(
  paste0("./", out_dir_top, "S01_ST_predictions/"),
  paste0("./", out_dir_top, "S02_ST_predictions/"),
  paste0("./", out_dir_top, "S89_ST_predictions/")
))
misty_res_lme3 <- collect_results(c(
  paste0("./", out_dir_top, "S04_ST_predictions/"),
  paste0("./", out_dir_top, "S06_ST_predictions/"),
  paste0("./", out_dir_top, "S87_ST_predictions/")
))

misty_list <- list(
  LME_1 = misty_res_lme1,
  LME_2 = misty_res_lme2,
  LME_3 = misty_res_lme3
)

spot_import_df <- purrr::map_dfr(names(misty_list), function(lme_name) {
  misty_list[[lme_name]]$importances.aggregated %>%
    dplyr::filter(
      Predictor %in% c("B.Mal", "CD4T.CXCL13", "FDC.CR2", "Mac.SPARC"),
      Target    %in% c("CD40LG", "TNFSF13B"),
      view      %in% c("intra", "juxta_2")
    ) %>%
    dplyr::mutate(subtype = lme_name)
})

spot_import_df$pred_subtype <- paste(spot_import_df$Predictor, spot_import_df$subtype, sep = "_")
max_abs_top <- max(spot_import_df$Importance, na.rm = TRUE)

# vline positions
x_interval   <- length(unique(spot_import_df$subtype))
x_intercepts <- seq(
  0.5 + x_interval,
  length(unique(spot_import_df$pred_subtype)) - 0.5,
  by = x_interval
)

figS7d_top_dotplot <- ggplot(spot_import_df, aes(x = pred_subtype, y = Target)) +
  geom_point(aes(size = Importance, color = Importance)) +
  facet_wrap(~view) +
  scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, max_abs_top)) +
  scale_size(range = c(0.5, 2), limits = c(0, max_abs_top)) +
  geom_vline(xintercept = x_intercepts, linetype = "dashed", color = "grey60", linewidth = 0.116) +
  theme_bw(base_size = 5) +
  labs(x = NULL, y = NULL) +
  theme(
    text            = element_text(size = 5),
    panel.grid.major= element_line(colour = "grey90", linewidth = 0.116),
    panel.grid.minor= element_blank(),
    panel.background= element_blank(),
    panel.border    = element_rect(colour = "black", linewidth = 0.116),
    axis.line       = element_line(colour = "black", linewidth = 0.116),
    axis.ticks      = element_line(colour = "black", linewidth = 0.116),
    axis.text.x     = element_text(colour = "black", angle = 60, vjust = 1, hjust = 1, size = 5),
    axis.text.y     = element_text(colour = "black", size = 5),
    strip.background= element_blank(),
    legend.position = "right",
    legend.title    = element_text(size = 5)
  )

ggsave(
  "figS7d_top.pdf", width = 12, height = 10, units = "cm",
  egg::set_panel_size(
    figS7d_top_dotplot,
    width  = unit(length(unique(spot_import_df$Predictor)) * 3 * 0.33, "cm"),
    height = unit(length(unique(spot_import_df$Target)) * 0.33, "cm")
  )
)

# ---- figS7d bottom (pathway-cell abundance) ----
out_dir_bottom <- "figS7d_bottom_"

for (slide_id in slide_ids) {
  slide <- subset(seurat_spatial, orig.ident == slide_id)
  
  cell_mtx    <- t(cell2loc_df[colnames(slide), ])
  progeny_mat <- progeny(
    slide, scale = FALSE, organism = "Human", top = 500, perm = 1,
    assay_name = "Spatial", return_assay = FALSE
  )
  pathway_mtx <- t(progeny_mat)
  merge_mtx   <- rbind(cell_mtx, pathway_mtx)
  
  slide@assays$predictions <- CreateAssayObject(counts = merge_mtx)
  DefaultAssay(slide)      <- "predictions"
  slide@images             <- seurat_spatial@images[slide_id]
  
  useful_features <- rownames(slide)
  
  invisible(run_colocalization(
    slide            = slide,
    useful_features  = useful_features,
    out_label        = slide_id,
    assay            = "predictions",
    misty_out_alias  = out_dir_bottom
  ))
}

misty_res_lme1 <- collect_results(c(
  paste0("./", out_dir_bottom, "S31_ST_predictions/"),
  paste0("./", out_dir_bottom, "S32_ST_predictions/"),
  paste0("./", out_dir_bottom, "S71_ST_predictions/")
))
misty_res_lme2 <- collect_results(c(
  paste0("./", out_dir_bottom, "S01_ST_predictions/"),
  paste0("./", out_dir_bottom, "S02_ST_predictions/"),
  paste0("./", out_dir_bottom, "S89_ST_predictions/")
))
misty_res_lme3 <- collect_results(c(
  paste0("./", out_dir_bottom, "S04_ST_predictions/"),
  paste0("./", out_dir_bottom, "S06_ST_predictions/"),
  paste0("./", out_dir_bottom, "S87_ST_predictions/")
))

misty_list <- list(
  LME_1 = misty_res_lme1,
  LME_2 = misty_res_lme2,
  LME_3 = misty_res_lme3
)

spot_import_df <- purrr::map_dfr(names(misty_list), function(lme_name) {
  misty_list[[lme_name]]$importances.aggregated %>%
    dplyr::filter(
      Predictor %in% c("Neu.G0S2", "Mac.MARCO"),
      Target    %in% c("Endo.ACKR1", "Endo.CXCL12", "VEGF"),
      view      %in% c("intra", "juxta_2")
    ) %>%
    dplyr::mutate(subtype = lme_name)
})

spot_import_df$pred_subtype <- paste(spot_import_df$Predictor, spot_import_df$subtype, sep = "_")
max_abs_bottom <- max(spot_import_df$Importance, na.rm = TRUE)

x_interval   <- length(unique(spot_import_df$subtype))
x_intercepts <- seq(
  0.5 + x_interval,
  length(unique(spot_import_df$pred_subtype)) - 0.5,
  by = x_interval
)

figS7d_bottom_dotplot <- ggplot(spot_import_df, aes(x = pred_subtype, y = Target)) +
  geom_point(aes(size = Importance, color = Importance)) +
  facet_wrap(~view) +
  scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, max_abs_bottom)) +
  scale_size(range = c(0.5, 2), limits = c(0, max_abs_bottom)) +
  geom_vline(xintercept = x_intercepts, linetype = "dashed", color = "grey60", linewidth = 0.116) +
  theme_bw(base_size = 5) +
  labs(x = NULL, y = NULL) +
  theme(
    text            = element_text(size = 5),
    panel.grid.major= element_line(colour = "grey90", linewidth = 0.116),
    panel.grid.minor= element_blank(),
    panel.background= element_blank(),
    panel.border    = element_rect(colour = "black", linewidth = 0.116),
    axis.line       = element_line(colour = "black", linewidth = 0.116),
    axis.ticks      = element_line(colour = "black", linewidth = 0.116),
    axis.text.x     = element_text(colour = "black", angle = 60, vjust = 1, hjust = 1, size = 5),
    axis.text.y     = element_text(colour = "black", size = 5),
    strip.background= element_blank(),
    legend.position = "right",
    legend.title    = element_text(size = 5)
  )

ggsave(
  "figS7d_bottom.pdf", width = 10, height = 10, units = "cm",
  egg::set_panel_size(
    figS7d_bottom_dotplot,
    width  = unit(length(unique(spot_import_df$Predictor)) * 3 * 0.33, "cm"),
    height = unit(length(unique(spot_import_df$Target)) * 0.33, "cm")
  )
)

# ---- figS7b ----
out_dir_S7b <- "figS7b_"

for (slide_id in slide_ids) {
  slide <- subset(seurat_spatial, orig.ident == slide_id)
  
  cell_mtx    <- t(cell2loc_df[colnames(slide), ])
  
  slide@assays$predictions <- CreateAssayObject(counts = cell_mtx)
  DefaultAssay(slide)      <- "predictions"
  slide@images             <- seurat_spatial@images[slide_id]
  
  useful_features <- rownames(slide)
  
  invisible(run_colocalization(
    slide            = slide,
    useful_features  = useful_features,
    out_label        = slide_id,
    assay            = "predictions",
    misty_out_alias  = out_dir_S7b
  ))
}

result_filedir <- grep(out_dir_S7b, list.dirs(), value = T)

misty_res <- collect_results(result_filedir)

pdf("figS7b_left.pdf")
mistyR::plot_interaction_communities(misty_res, "intra", cutoff = 0.5)
dev.off()

pdf("figS7b_right.pdf")
mistyR::plot_interaction_communities(misty_res, "juxta_2", cutoff = 0.5)
dev.off()
