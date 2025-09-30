# ---- packages ----
library(Seurat)
library(AUCell)
library(tidyverse)
library(GSVA)
library(GSEABase)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(egg)
library(grid)

# ---- paths & sources ----
snrna_data_dir <- "../../data/snrna_rnaseq"
bulk_data_dir  <- "../../data/bulk_rnaseq"
src_dir        <- "../source"
fig_dir        <- "."

source(file.path(src_dir, "custom_colors.R"))

# ---- snRNA AUCell: inputs ----
seurat_b    <- readRDS(file.path(snrna_data_dir, "snrna_b_subset_seurat.rds"))
seurat_bmal <- subset(seurat_b, fine_cell_type == "B_Mal")
expr_counts <- as.matrix(seurat_bmal@assays$RNA@counts)

# ---- snRNA AUCell: rankings & AUC ----
cell_rankings  <- AUCell_buildRankings(expr_counts, nCores = 12, plotStats = TRUE)

gene_sets_df   <- read.csv(file.path(src_dir, "MP_list.csv"),
                           header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gene_sets_list <- as.list(gene_sets_df)

cells_auc     <- AUCell_calcAUC(gene_sets_list, cell_rankings)
saveRDS(cells_auc, file.path(fig_dir, "cells_AUC.rds"))

auc_score_mat <- AUCell::getAUC(cells_auc)  # gene sets x cells
stopifnot(identical(rownames(seurat_bmal@meta.data), colnames(auc_score_mat)))
seurat_bmal@meta.data <- cbind(seurat_bmal@meta.data, t(auc_score_mat))

# ---- snRNA AUCell: per-MP boxplots ----
boxplot_input <- seurat_bmal@meta.data %>%
  dplyr::select(subtype, tidyselect::matches("^MP_\\d+$")) %>%
  tidyr::pivot_longer(cols = -subtype, names_to = "MP", values_to = "score")

get_boxplot_limits <- function(df) {
  p  <- ggplot(df, aes(x = subtype, y = score)) + geom_boxplot()
  s  <- ggplot_build(p)$data[[1]]
  tibble::tibble(ymin = min(s$ymin, na.rm = TRUE), ymax = max(s$ymax, na.rm = TRUE))
}


limits_by_mp <- boxplot_input %>%
  split(.$MP) %>%
  purrr::imap_dfr(function(df, nm) {
    lim <- get_boxplot_limits(df)
    dplyr::mutate(lim, MP = nm, .before = 1)
  })

snrna_mp_plots <- boxplot_input %>%
  split(.$MP) %>%
  lapply(function(df_mp) {
    y_range <- limits_by_mp %>% filter(MP == unique(df_mp$MP))
    
    ggplot(df_mp, aes(x = subtype, y = score, fill = subtype)) +
      geom_boxplot(outlier.shape = NA, size = 0.116, width = 0.6) +
      scale_fill_manual(values = subtype_col1) +
      coord_cartesian(ylim = c(y_range$ymin, y_range$ymax)) +
      theme_minimal(base_size = 5) +
      theme(
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.line        = element_line(colour = "black", linewidth = 0.116),
        axis.ticks       = element_line(colour = "black", linewidth = 0.116),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        axis.text.y      = element_text(colour = "black", size = 5),
        legend.position  = "none",
        plot.title       = element_text(size = 5, hjust = 0.5)
      ) +
      ggtitle(unique(df_mp$MP)) +
      labs(x = NULL, y = NULL)
  })

snrna_mp_panel <- wrap_plots(snrna_mp_plots, nrow = 1)
ggsave(
  filename = file.path(fig_dir, "fig3h_top.pdf"),
  plot = snrna_mp_panel,
  width = 12, height = 2, units = "cm"
)

# ---- bulk GSVA: inputs ----
bulk_expr_path     <- file.path(bulk_data_dir, "bulk_expression_matrix.txt")
bulk_clinical_path <- file.path(bulk_data_dir, "bulk_clinical_metadata.txt")

clinical_df <- read.table(
  bulk_clinical_path, header = TRUE, sep = "\t", row.names = 1,
  check.names = FALSE, na.strings = c("", NA)
)
clinical_tumor_df <- clinical_df[clinical_df$Group == "Tumor", , drop = FALSE]

# ---- bulk GSVA: expression & scores ----
bulk_expr_mat <- read.table(bulk_expr_path, header = TRUE, row.names = 1, check.names = FALSE)
bulk_expr_mat <- bulk_expr_mat[, rownames(clinical_tumor_df), drop = FALSE]
bulk_expr_mat <- as.matrix(log2(bulk_expr_mat + 1))

gsva_score_mat <- gsva(bulk_expr_mat, gene_sets_list, method = "gsva", kcdf = "Gaussian")

# ---- bulk GSVA: clinical groupings ----
clinical_tumor_df <- clinical_tumor_df %>%
  mutate(
    Age_group   = dplyr::case_when(Age < 60 ~ "<60", Age >= 60 ~ ">=60"),
    Stage_group = dplyr::case_when(
      Ann_Arbor_stage %in% c("I", "II", "III") ~ "Stage I_III",
      Ann_Arbor_stage %in% c("IV")            ~ "Stage IV"
    ),
    Site_group  = dplyr::case_when(!Site %in% c("Ocular_adnexa", "Stomach") ~ "Others", TRUE ~ Site)
  )

# ---- bulk GSVA: helper to plot per-MP boxplots by group ----
plot_mp_box <- function(score_mat,
                        group_vec,
                        group_name = "group",
                        fill_cols = NULL,
                        out_file = NULL,
                        out_w = 16, out_h = 2, out_unit = "cm") {
  
  # 1) data wrangling
  long_df <- score_mat %>%
    t() %>%
    as.data.frame() %>%
    mutate(!!group_name := as.vector(group_vec)) %>%
    drop_na(!!sym(group_name)) %>%
    pivot_longer(
      cols = -all_of(group_name),
      names_to = "MP",
      values_to = "score"
    )
  
  # set factor order for groups (if color palette is provided)
  if (!is.null(fill_cols)) {
    long_df[[group_name]] <- factor(long_df[[group_name]], levels = names(fill_cols))
  } else {
    long_df[[group_name]] <- factor(long_df[[group_name]])
  }
  
  # 2) plot per MP
  plot_list <- long_df %>%
    group_split(MP) %>%
    map(function(df_mp) {
      ggplot(df_mp, aes(x = .data[[group_name]],
                        y = score,
                        fill = .data[[group_name]])) +
        geom_boxplot(outlier.shape = NA, size = 0.116, width = 0.6) +
        theme_minimal(base_size = 5) +
        theme(
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black", size = 0.116),
          axis.ticks = element_line(colour = "black", size = 0.116),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(colour = "black", size = 5),
          legend.position = "none",
          plot.title = element_text(size = 5, hjust = 0.5)
        ) +
        ggtitle(unique(df_mp$MP)) +
        coord_cartesian(clip = "off") +
        xlab("") + ylab("") +
        {
          if (!is.null(fill_cols)) scale_fill_manual(values = fill_cols)
        } +
        stat_pwc(
          aes(group = .data[[group_name]]),
          label.size = 5 * 25.4 / 72,
          tip.length = 0,
          size = 0.116,
          hide.ns = TRUE,
          vjust = 0.5,
          method = "wilcox.test",
          p.adjust.method = "BH",
          label = "p.adj.signif"
        )
    })
  
  # 3) assemble
  panel_wrap <- cowplot::plot_grid(plotlist = plot_list,
                                   nrow = 1,
                                   rel_widths = rep(1.5, length(plot_list)),
                                   align = "hv")
  
  legend <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = "right"))
  panel_wrap <- cowplot::plot_grid(panel_wrap, legend, ncol = 2, rel_widths = c(7, 1))
  
  # 4) optional save
  if (!is.null(out_file)) {
    ggsave(filename = out_file, plot = panel_wrap,
           width = out_w, height = out_h, units = out_unit)
  }
  
  return(panel_wrap)
}

# ---- snRNA-seq AUCell: outputs with significance ----
plot_mp_box(
  score_mat = auc_score_mat,
  group_vec = seurat_bmal@meta.data[, "subtype"],
  group_name = "subtype",
  fill_cols = subtype_col1,
  out_file = file.path(fig_dir, "fig3h_top_sig.pdf")
)

# ---- bulk GSVA: outputs ----
plot_mp_box(
  score_mat = gsva_score_mat,
  group_vec = clinical_tumor_df[, "LME_subtype"],
  group_name = "subtype",
  fill_cols = subtype_col1,
  out_file = file.path(fig_dir, "fig3h_bottom.pdf")
)

plot_mp_box(
  score_mat = gsva_score_mat,
  group_vec = clinical_tumor_df[, "Stage_group"],
  group_name = "stage",
  fill_cols = stage_group_col,
  out_file = file.path(fig_dir, "figS3b_stage.pdf")
)

plot_mp_box(
  score_mat = gsva_score_mat,
  group_vec = clinical_tumor_df[, "Age_group"],
  group_name = "age",
  fill_cols = age_group_col,
  out_file = file.path(fig_dir, "figS3b_age.pdf")
)

plot_mp_box(
  score_mat = gsva_score_mat,
  group_vec = clinical_tumor_df[, "BIRC3_MALT_fusion"],
  group_name = "fusion",
  fill_cols = fusion_col,
  out_file = file.path(fig_dir, "figS3b_fusion.pdf")
)

plot_mp_box(
  score_mat = gsva_score_mat,
  group_vec = clinical_tumor_df[, "Site_group"],
  group_name = "site",
  fill_cols = site_col,
  out_file = file.path(fig_dir, "figS3b_site.pdf")
)