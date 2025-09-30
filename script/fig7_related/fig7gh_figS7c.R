# ---- packages ----
library(Seurat)
library(ggpubr)
library(ggridges)
library(dplyr)
library(png)
library(jsonlite)

# ---- paths & sources ----
data_dir <- "../../data/spatial_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))
source(file.path(src_dir, "calculate_spatial_distance.R"))

# ---- discover samples ----
sample_ids <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)

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

# ---- per-slide pairwise distance comparison ----
global_seed <- 421
slide_ids   <- sample_ids

slide_result_list <- list()

for (slide_id in slide_ids) {
  seu_slide <- subset(seurat_spatial, orig.ident == slide_id)
  
  # microns-per-pixel conversion using fiducial diameter reported by 10x
  fiducial_diam <- seurat_spatial@images[[slide_id]]@scale.factors$fiducial
  
  # positions in pixels -> convert to microns (65 µm nominal spot diameter on Visium)
  pos_df <- seurat_spatial@images[[slide_id]]@coordinates
  pos_df$imagerow <- pos_df$imagerow * (65 / fiducial_diam)
  pos_df$imagecol <- pos_df$imagecol * (65 / fiducial_diam)
  
  coords_mat <- as.matrix(pos_df[, c("imagerow", "imagecol")])
  dist_mat   <- as.matrix(dist(coords_mat))
  
  prop_mat <- cell2loc_df[rownames(pos_df), , drop = FALSE]
  celltypes <- colnames(prop_mat)
  
  pair_result <- list()
  
  for (i in seq_along(celltypes)) {
    for (j in seq_along(celltypes)) {
      ct1 <- celltypes[i]
      ct2 <- celltypes[j]
      if (ct1 == ct2) next
      
      # two-means threshold per cell type (reproducible with fixed seed)
      set.seed(global_seed)
      km1 <- kmeans(prop_mat[[ct1]], centers = 2, nstart = 10)
      t1  <- mean(range(tapply(prop_mat[[ct1]], km1$cluster, mean)))
      
      set.seed(global_seed)
      km2 <- kmeans(prop_mat[[ct2]], centers = 2, nstart = 10)
      t2  <- mean(range(tapply(prop_mat[[ct2]], km2$cluster, mean)))
      
      # distance comparison (function defined in calculate_spatial_distance.R)
      cmp_df <- distance_comparison(
        distances   = dist_mat,
        proportions = prop_mat,
        topic1      = ct1,
        topic2      = ct2,
        t1_thresh   = t1,
        t2_thresh   = t2
      )
      
      pair_name <- paste0(ct1, "_vs_", ct2)
      pair_result[[pair_name]] <- cmp_df
    }
  }
  
  slide_result_list[[slide_id]] <- do.call(rbind, pair_result)
}

distance_df <- do.call(rbind, slide_result_list)

# add slide id (assumes rownames formatted as "<slide>.<rowindex>")
distance_df$slide_id <- sapply(strsplit(rownames(distance_df), "[.]"), `[`, 1)

# annotate LME subtype by slide id
distance_df <- mutate(
  distance_df,
  subtype = case_when(
    slide_id %in% c("S31_ST", "S32_ST", "S71_ST") ~ "LME_1",
    slide_id %in% c("S01_ST", "S02_ST", "S89_ST") ~ "LME_2",
    slide_id %in% c("S04_ST", "S06_ST", "S87_ST") ~ "LME_3",
    TRUE ~ NA_character_
  )
)

# ---- Fig 7g: B cell niche distances by subtype ----
plot_df <- subset(
  distance_df,
  barcodes %in% c(
    "B_Mal vs. CD8T_HAVCR2",
    "B_Mal vs. CD4T_CXCL13",
    "B_Mal vs. FDC_CR2",
    "B_Mal vs. Mac_SPARC"
  )
)
plot_df$barcodes <- factor(plot_df$barcodes)
plot_df$subtype  <- factor(plot_df$subtype)

median_df <- plot_df |>
  group_by(barcodes, subtype) |>
  summarise(med = median(min), .groups = "drop")

ridge_p <- ggplot(plot_df, aes(x = min / 1000, y = subtype, fill = subtype)) +
  geom_density_ridges(alpha = 0.5, scale = 1, size = 0.116) +
  geom_vline(
    data = median_df,
    aes(xintercept = med / 1000, colour = subtype),
    linetype = "dashed", linewidth = 0.116, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = subtype_col1) +
  scale_color_manual(values = subtype_col1) +
  labs(x = "Distance (mm)", y = NULL) +
  facet_wrap(~ barcodes, ncol = 4, scales = "free_x", switch = "y") +
  theme_bw() +
  theme(
    text            = element_text(size = 5),
    axis.text.x     = element_text(colour = "black", size = 5),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.line       = element_line(colour = "black", linewidth = 0.116),
    axis.ticks.x    = element_line(colour = "black", linewidth = 0.116),
    legend.title    = element_text(size = 5),
    strip.background= element_blank(),
    strip.text      = element_text(size = 5),
    panel.grid.major= element_blank(),
    panel.grid.minor= element_blank(),
    panel.border    = element_rect(colour = "black", linewidth = 0.116)
  )

ggsave(
  filename = file.path(fig_dir, "fig7g.pdf"),
  plot     = egg::set_panel_size(ridge_p, width = unit(1, "cm"), height = unit(1, "cm")),
  width = 10, height = 10, units = "cm"
)

# ---- Fig S7c (left): Endo_CXCL12 vs myeloid subsets ----
plot_df <- subset(
  distance_df,
  barcodes %in% c(
    "Endo_CXCL12 vs. Neu_G0S2",
    "Endo_CXCL12 vs. Mac_MARCO",
    "Endo_CXCL12 vs. Mac_SPARC",
    "Endo_CXCL12 vs. Mac_CHIT1",
    "Endo_CXCL12 vs. Mac_STAB1"
  )
)
plot_df$cell_type <- sapply(strsplit(plot_df$barcodes, " "), `[`, 3)
plot_df$cell_type <- factor(
  plot_df$cell_type,
  levels = c("Neu_G0S2", "Mac_MARCO", "Mac_SPARC", "Mac_STAB1", "Mac_CHIT1")
)
plot_df$subtype <- factor(plot_df$subtype)

median_df <- plot_df |>
  group_by(cell_type) |>
  summarise(med = median(min), .groups = "drop")

ridge_p <- ggplot(plot_df, aes(x = min / 1000, y = cell_type, fill = cell_type)) +
  geom_density_ridges(alpha = 0.5, scale = 1, size = 0.116) +
  geom_vline(
    data = median_df,
    aes(xintercept = med / 1000, colour = cell_type),
    linetype = "dashed", linewidth = 0.116
  ) +
  scale_fill_manual(values = fine_cell_type_col) +
  scale_color_manual(values = fine_cell_type_col) +
  labs(x = "Distance (mm)", y = NULL) +
  theme_bw() +
  theme(
    text            = element_text(size = 5),
    axis.text.x     = element_text(colour = "black", size = 5),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.line       = element_line(colour = "black", linewidth = 0.116),
    axis.ticks.x    = element_line(colour = "black", linewidth = 0.116),
    legend.title    = element_text(size = 5),
    strip.background= element_blank(),
    strip.text      = element_text(size = 5),
    panel.grid.major= element_blank(),
    panel.grid.minor= element_blank(),
    panel.border    = element_rect(colour = "black", linewidth = 0.116)
  ) +
  scale_y_discrete(limits = rev)

ggsave(
  filename = file.path(fig_dir, "figS7c_left.pdf"),
  plot     = egg::set_panel_size(ridge_p, width = unit(1, "cm"), height = unit(2, "cm")),
  width = 10, height = 10, units = "cm"
)

# ---- Fig S7c (right): Endo_ACKR1 vs myeloid subsets ----
plot_df <- subset(
  distance_df,
  barcodes %in% c(
    "Endo_ACKR1 vs. Neu_G0S2",
    "Endo_ACKR1 vs. Mac_MARCO",
    "Endo_ACKR1 vs. Mac_SPARC",
    "Endo_ACKR1 vs. Mac_CHIT1",
    "Endo_ACKR1 vs. Mac_STAB1"
  )
)
plot_df$cell_type <- sapply(strsplit(plot_df$barcodes, " "), `[`, 3)
plot_df$cell_type <- factor(
  plot_df$cell_type,
  levels = c("Neu_G0S2", "Mac_MARCO", "Mac_SPARC", "Mac_STAB1", "Mac_CHIT1")
)
plot_df$subtype <- factor(plot_df$subtype)

median_df <- plot_df |>
  group_by(cell_type) |>
  summarise(med = median(min), .groups = "drop")

ridge_p <- ggplot(plot_df, aes(x = min / 1000, y = cell_type, fill = cell_type)) +
  geom_density_ridges(alpha = 0.5, scale = 1, size = 0.116) +
  geom_vline(
    data = median_df,
    aes(xintercept = med / 1000, colour = cell_type),
    linetype = "dashed", linewidth = 0.116
  ) +
  scale_fill_manual(values = fine_cell_type_col) +
  scale_color_manual(values = fine_cell_type_col) +
  labs(x = "Distance (mm)", y = NULL) +
  theme_bw() +
  theme(
    text            = element_text(size = 5),
    axis.text.x     = element_text(colour = "black", size = 5),
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    axis.line       = element_line(colour = "black", linewidth = 0.116),
    axis.ticks.x    = element_line(colour = "black", linewidth = 0.116),
    legend.title    = element_text(size = 5),
    strip.background= element_blank(),
    strip.text      = element_text(size = 5),
    panel.grid.major= element_blank(),
    panel.grid.minor= element_blank(),
    panel.border    = element_rect(colour = "black", linewidth = 0.116)
  ) +
  scale_y_discrete(limits = rev)

ggsave(
  filename = file.path(fig_dir, "figS7c_right.pdf"),
  plot     = egg::set_panel_size(ridge_p, width = unit(1, "cm"), height = unit(2, "cm")),
  width = 10, height = 10, units = "cm"
)

# ---- per-slide minimum distance to a single topic (cell type) ----
topic_results_by_slide <- list()

for (sid in slide_ids) {
  seu_slide <- subset(seurat_spatial, orig.ident == sid)
  
  # pixel -> micron conversion
  fiducial_diam <- seurat_spatial@images[[sid]]@scale.factors$fiducial
  pos_df <- seurat_spatial@images[[sid]]@coordinates
  pos_df$imagerow <- pos_df$imagerow * (65 / fiducial_diam)
  pos_df$imagecol <- pos_df$imagecol * (65 / fiducial_diam)
  
  coords_mat <- as.matrix(pos_df[, c("imagerow", "imagecol")])
  dist_mat   <- as.matrix(dist(coords_mat))
  
  prop_mat   <- cell2loc_df[rownames(pos_df), , drop = FALSE]
  celltypes  <- colnames(prop_mat)
  
  topic_result_list <- list()
  
  for (ct in celltypes) {
    set.seed(global_seed)
    km   <- kmeans(prop_mat[[ct]], centers = 2, nstart = 10)
    thr  <- mean(range(tapply(prop_mat[[ct]], km$cluster, mean)))
    
    res_df <- min_distance_to_topic(
      distances   = dist_mat,
      proportions = prop_mat,
      topic       = ct,
      thresh      = thr
    )
    
    topic_result_list[[ct]] <- res_df
  }
  
  topic_results_by_slide[[sid]] <- do.call(rbind, topic_result_list)
}

distance_topic_df <- do.call(rbind, topic_results_by_slide)

# annotate slide + subtype
distance_topic_df$slide_id <- sapply(strsplit(rownames(distance_topic_df), "[.]"), `[`, 1)
distance_topic_df <- dplyr::mutate(
  distance_topic_df,
  subtype = dplyr::case_when(
    slide_id %in% c("S31_ST", "S32_ST", "S71_ST") ~ "LME_1",
    slide_id %in% c("S01_ST", "S02_ST", "S89_ST") ~ "LME_2",
    slide_id %in% c("S04_ST", "S06_ST", "S87_ST") ~ "LME_3",
    TRUE ~ NA_character_
  )
)

# ---- gene–distance correlation (per topic & gene) ----
target_topics <- c(
  "CD4T_CXCL13", "Mac_SPARC", "FDC_CR2",
  "Neu_G0S2", "Mac_MARCO", "Endo_CXCL12", "Endo_ACKR1"
)
target_genes <- c(
  "CD40LG", "VEGFA"
)

gene_distance_rows <- list()

for (sid in unique(distance_topic_df$slide_id)) {
  seu_slide <- subset(seurat_spatial, orig.ident == sid)
  seu_slide <- NormalizeData(seu_slide, verbose = FALSE)
  
  for (topic in target_topics) {
    dsub <- subset(distance_topic_df, slide_id == sid & reference_topic == topic)
    barcodes <- dsub$spot
    barcodes <- barcodes[barcodes %in% colnames(seu_slide)]
    if (length(barcodes) == 0) next
    
    for (g in target_genes) {
      if (!(g %in% rownames(seu_slide@assays$Spatial@data))) {
        warning(paste("Gene", g, "not found in slide", sid))
        next
      }
      expr_vec <- seu_slide@assays$Spatial@data[g, barcodes]
      dist_vec <- dsub[match(barcodes, dsub$spot), "min_distance_to_topic"]
      
      gene_distance_rows[[length(gene_distance_rows) + 1]] <- data.frame(
        slide_id        = sid,
        reference_topic = topic,
        gene            = g,
        spot            = barcodes,
        expr            = as.numeric(expr_vec),
        distance        = as.numeric(dist_vec),
        stringsAsFactors = FALSE
      )
    }
  }
}

gene_distance_df <- do.call(rbind, gene_distance_rows)

gene_distance_df <- dplyr::mutate(
  gene_distance_df,
  subtype = dplyr::case_when(
    slide_id %in% c("S31_ST", "S32_ST", "S71_ST") ~ "LME_1",
    slide_id %in% c("S01_ST", "S02_ST", "S89_ST") ~ "LME_2",
    slide_id %in% c("S04_ST", "S06_ST", "S87_ST") ~ "LME_3",
    TRUE ~ NA_character_
  )
)

cor_summary_df <- gene_distance_df %>%
  dplyr::group_by(subtype, reference_topic, gene) %>%
  dplyr::summarise(
    spearman_rho = cor(expr, distance, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(gene_topic = paste(gene, reference_topic, sep = " | "))

# ---- Fig 7h (CD40LG) ----
cor_plot_df <- cor_summary_df
cor_plot_df$gene_topic <- factor(cor_plot_df$gene_topic,
                                 levels = c("CD40LG | CD4T_CXCL13", "CD40LG | Mac_SPARC", "CD40LG | FDC_CR2",     
                                            "VEGFA | Neu_G0S2", "VEGFA | Mac_MARCO", "VEGFA | Endo_ACKR1", "VEGFA | Endo_CXCL12"))

cor_plot_df_sub <- subset(
  cor_plot_df,
  reference_topic %in% c("CD4T_CXCL13", "Mac_SPARC", "FDC_CR2") & gene == "CD40LG"
)

max_abs_rho <- max(abs(cor_plot_df$spearman_rho), na.rm = TRUE)

cd40_dot <- ggplot(
  cor_plot_df_sub,
  aes(x = gene_topic, y = subtype, size = abs(spearman_rho), color = spearman_rho)
) +
  geom_point() +
  scale_color_gradientn(colors = rev(heatmap_col1),
                        limits = c(-max_abs_rho, max_abs_rho)) +
  scale_size(range = c(0, 2), limits = c(0, max_abs_rho)) +
  labs(x = "", y = "", size = "|rho|", color = "Spearman rho") +
  theme_bw(base_size = 5) +
  theme(
    text = element_text(size = 5),
    panel.grid.major = element_line(colour = "grey90", size = 0.116),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", size = 0.116),
    axis.line = element_line(colour = "black", size = 0.116),
    axis.ticks = element_line(colour = "black", size = 0.116),
    axis.text.x = element_text(colour = "black", angle = 30, vjust = 1, hjust = 1, size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    legend.position = "right",
    legend.title = element_text(size = 5)
  ) +
  guides(
    size  = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(stroke = 0.4)),
    color = guide_colourbar(title.position = "top", title.hjust = 0.5)
  )

ggsave(
  filename = file.path(fig_dir, "fig7h_CD40LG.pdf"),
  plot     = egg::set_panel_size(
    cd40_dot,
    width  = unit(length(unique(cor_plot_df_sub$gene_topic)) * 0.33, "cm"),
    height = unit(length(unique(cor_plot_df_sub$subtype)) * 0.33, "cm")
  ),
  width = 10, height = 10, units = "cm"
)

# ---- Fig 7h (VEGF/endothelium–myeloid axis topics) ----
cor_plot_df_sub <- subset(
  cor_plot_df,
  reference_topic %in% c("Neu_G0S2", "Mac_MARCO", "Endo_CXCL12", "Endo_ACKR1") & gene == "VEGFA"
)

vegf_dot <- ggplot(
  cor_plot_df_sub,
  aes(x = gene_topic, y = subtype, size = abs(spearman_rho), color = spearman_rho)
) +
  geom_point() +
  scale_color_gradientn(colors = rev(heatmap_col1),
                        limits = c(-max_abs_rho, max_abs_rho)) +
  scale_size(range = c(0, 2), limits = c(0, max_abs_rho)) +
  labs(x = "", y = "", size = "|rho|", color = "Spearman rho") +
  theme_bw(base_size = 5) +
  theme(
    text = element_text(size = 5),
    panel.grid.major = element_line(colour = "grey90", size = 0.116),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", size = 0.116),
    axis.line = element_line(colour = "black", size = 0.116),
    axis.ticks = element_line(colour = "black", size = 0.116),
    axis.text.x = element_text(colour = "black", angle = 30, vjust = 1, hjust = 1, size = 5),
    axis.text.y = element_text(colour = "black", size = 5),
    legend.position = "right",
    legend.title = element_text(size = 5)
  ) +
  guides(
    size  = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(stroke = 0.4)),
    color = guide_colourbar(title.position = "top", title.hjust = 0.5)
  )

ggsave(
  filename = file.path(fig_dir, "fig7h_VEGFA.pdf"),
  plot     = egg::set_panel_size(
    vegf_dot,
    width  = unit(length(unique(cor_plot_df_sub$gene_topic)) * 0.33, "cm"),
    height = unit(length(unique(cor_plot_df_sub$subtype)) * 0.33, "cm")
  ),
  width = 10, height = 10, units = "cm"
)