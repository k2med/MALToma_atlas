# ---- packages ----
library(Seurat)
library(reshape2)
library(NMF)
library(ggplot2)
library(scales)
library(tidyverse)

# ---- paths & sources ----
data_dir <- "../../data/snrna_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

source(file.path(src_dir, "custom_colors.R"))
source(file.path(src_dir, "robust_nmf_programs.R"))

# ---- load & subset ----
b_cell_seurat  <- readRDS(file.path(data_dir, "snrna_b_subset_seurat.rds"))
b_mal_seurat   <- subset(b_cell_seurat, fine_cell_type == "B_Mal")

samples_use <- unique(b_mal_seurat$orig.ident)

# ---- per-sample NMF (collect W) ----
genes_nmf_w_basis <- lapply(samples_use, function(sid) {
  b_mal_seurat_sub <- b_mal_seurat %>%
    subset(orig.ident == sid) %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData(vars.to.regress = c("percent.mt", "nCount_RNA", "nFeature_RNA"))
  
  exp_sub <- b_mal_seurat_sub@assays$RNA@scale.data
  exp_sub[exp_sub < 0] <- 0
  
  nmf_sub <- nmf(x = exp_sub, rank = 4:9, method = "snmf/r", nrun = 10)
  
  genes_nmf_w_sub_list <- lapply(4:9, function(rk) {
    w_mat <- nmf_sub$fit[[as.character(rk)]]@fit@W
    colnames(w_mat) <- paste0(sid, "_rank4_9_nruns10.", rk, ".", seq_len(rk))
    w_mat
  })
  
  genes_nmf_w_sub <- Reduce(function(m, n) cbind(m, n), genes_nmf_w_sub_list)
  
  saveRDS(genes_nmf_w_sub, paste0(sid, "_rank4_9_nruns10.RDS"))
  genes_nmf_w_sub
})

names(genes_nmf_w_basis) <- paste0(samples_use, "_rank4_9_nruns10")

# ---- select NMF programs ----
# parameters
intra_min_param <- 35
inter_min_param <- 10
intra_max_param <- 10

# top 50 genes per program
nmf_programs <- lapply(genes_nmf_w_basis, function(w)
  apply(w, 2, function(v) names(sort(v, decreasing = TRUE))[1:50]))
nmf_programs <- lapply(nmf_programs, toupper)

# robust program selection per sample; deduplicate across ranks; inter-sample filter
nmf_filter <- robust_nmf_programs(
  nmf_programs,
  intra_min   = intra_min_param,
  intra_max   = intra_max_param,
  inter_filter= TRUE,
  inter_min   = inter_min_param
)

nmf_programs <- lapply(nmf_programs, function(x) x[, colnames(x) %in% nmf_filter, drop = FALSE])
nmf_programs <- do.call(cbind, nmf_programs)

# similarity between programs (overlap sizes)
nmf_intersect <- apply(nmf_programs, 2, function(x)
  apply(nmf_programs, 2, function(y) length(intersect(x, y))))

# hierarchical clustering of similarity matrix
nmf_intersect_hc <- hclust(as.dist(50 - nmf_intersect), method = "average")
nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect    <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

# ---- cluster selected NMF programs to generate MPs ----
# clustering parameters
min_intersect_initial <- 10  # the minimal intersection cutoff for defining the first NMF program in a cluster
min_intersect_cluster <- 10  # the minimal intersection cutoff for adding a new NMF to the forming cluster
min_group_size        <- 3   # the minimal group size to consider for defining the first NMF program in a cluster

sorted_intersection <- sort(
  apply(nmf_intersect, 2, function(x) length(which(x >= min_intersect_initial)) - 1),
  decreasing = TRUE
)

cluster_list <- list()  # NMF ids per cluster
mp_list      <- list()  # 50-gene MP per cluster
k            <- 1
curr_cluster <- character(0)

nmf_intersect_original <- nmf_intersect

while (sorted_intersection[1] > min_group_size) {
  
  curr_cluster <- c(curr_cluster, names(sorted_intersection[1]))
  
  # seed genes for current MP
  genes_mp <- nmf_programs[, names(sorted_intersection[1])]
  nmf_programs <- nmf_programs[, -match(names(sorted_intersection[1]), colnames(nmf_programs))]
  inter_with_mp <- sort(apply(nmf_programs, 2, function(x) length(intersect(genes_mp, x))), decreasing = TRUE)
  nmf_history   <- genes_mp
  
  while (inter_with_mp[1] >= min_intersect_cluster) {
    
    curr_cluster <- c(curr_cluster, names(inter_with_mp)[1])
    
    genes_mp_tmp   <- sort(table(c(nmf_history, nmf_programs[, names(inter_with_mp)[1]])), decreasing = TRUE)
    genes_at_border <- genes_mp_tmp[which(genes_mp_tmp == genes_mp_tmp[50])]
    
    if (length(genes_at_border) > 1) {
      # break ties by NMF scores across programs in current cluster
      genes_curr_scores <- c()
      for (nid in curr_cluster) {
        curr_study <- paste(strsplit(nid, "[.]")[[1]][1], collapse = ".")
        Q <- genes_nmf_w_basis[[curr_study]][
          match(names(genes_at_border), toupper(rownames(genes_nmf_w_basis[[curr_study]])))[
            !is.na(match(names(genes_at_border), toupper(rownames(genes_nmf_w_basis[[curr_study]]))))
          ],
          nid
        ]
        names(Q) <- names(genes_at_border[!is.na(match(
          names(genes_at_border),
          toupper(rownames(genes_nmf_w_basis[[curr_study]]))
        ))])
        genes_curr_scores <- c(genes_curr_scores, Q)
      }
      genes_curr_scores <- genes_curr_scores[!is.na(genes_curr_scores)]
      genes_curr_scores_sort <- sort(genes_curr_scores, decreasing = TRUE)
      genes_curr_scores_sort <- genes_curr_scores_sort[unique(names(genes_curr_scores_sort))]
      
      genes_mp_tmp <- c(
        names(genes_mp_tmp[which(genes_mp_tmp > genes_mp_tmp[50])]),
        names(genes_curr_scores_sort)
      )
    } else {
      genes_mp_tmp <- names(genes_mp_tmp)[1:50]
    }
    
    nmf_history <- c(nmf_history, nmf_programs[, names(inter_with_mp)[1]])
    genes_mp    <- genes_mp_tmp[1:50]
    
    nmf_programs <- nmf_programs[, -match(names(inter_with_mp)[1], colnames(nmf_programs))]
    inter_with_mp <- sort(apply(nmf_programs, 2, function(x) length(intersect(genes_mp, x))), decreasing = TRUE)
  }
  
  cluster_list[[paste0("Cluster_", k)]] <- curr_cluster
  mp_list[[paste0("MP_", k)]]           <- genes_mp
  k <- k + 1
  
  # remove clustered programs and re-score remaining
  nmf_intersect <- nmf_intersect[
    -match(curr_cluster, rownames(nmf_intersect)),
    -match(curr_cluster, colnames(nmf_intersect))
  ]
  
  sorted_intersection <- sort(
    apply(nmf_intersect, 2, function(x) length(which(x >= min_intersect_initial)) - 1),
    decreasing = TRUE
  )
  
  curr_cluster <- character(0)
  print(dim(nmf_intersect)[2])
}

# ---- order similarity matrix by new clusters and plot ----
inds_sorted <- c()
for (j in seq_along(cluster_list)) {
  inds_sorted <- c(inds_sorted, match(cluster_list[[j]], colnames(nmf_intersect_original)))
}
inds_new <- c(inds_sorted, which(is.na(match(seq_len(ncol(nmf_intersect_original)), inds_sorted))))

# remove the MP suspected to reflect doublets
inds_new_sub <- inds_new[c(1:10, 20:length(inds_new))]

nmf_intersect_melt_new <- reshape2::melt(nmf_intersect_original[inds_new_sub, inds_new_sub])

heatmap_plot <- ggplot(nmf_intersect_melt_new, aes(x = Var1, y = Var2, fill = value, color = value)) +
  geom_tile() +
  scale_color_gradient2(
    limits = c(2, 25),
    low = heatmap_col4[1:111], mid = heatmap_col4[112:222], high = heatmap_col4[223:333],
    midpoint = 13.5, oob = squish, name = "Overlap size"
  ) +
  scale_fill_gradient2(
    limits = c(2, 25),
    low = heatmap_col4[1:111], mid = heatmap_col4[112:222], high = heatmap_col4[223:333],
    midpoint = 13.5, oob = squish, name = "Overlap size"
  ) +
  theme(
    axis.ticks = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.background = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.text.align = 0.5,
    legend.justification = "bottom",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  scale_y_discrete(limits=rev) +
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

ggsave(
  "fig3a.pdf",
  egg::set_panel_size(heatmap_plot, width = unit(4, "cm"), height = unit(4, "cm")),
  width = 9, height = 6, units = "cm"
)