# ---- running CytoSig in a Docker ----
docker run -it -w /tests -v "$(pwd):/tests" data2intelligence/data2intelligence-suite
CytoSig_run.py -i bulk_expression_matrix.txt -o output -e 1

# ---- packages ----
library(ggpubr)
library(ggplot2)
library(dplyr)
library(egg)
library(grid)

# ---- paths & sources ----
data_dir   <- "../../data/bulk_rnaseq"
src_dir    <- "../source"
fig_dir    <- "."

# files used in this step
bulk_clinical_path <- file.path(data_dir, "bulk_clinical_metadata.txt")
cytosig_out_path    <- "output.Zscore"

source(file.path(src_dir, "custom_colors.R"))

# ---- load CytoSig results ----
cytosig_df <- read.table(cytosig_out_path, header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE)

# ---- clinical (tumor only) ----
clinical_df <- read.table(
  bulk_clinical_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)
clinical_tumor_df <- clinical_df[clinical_df$Group == "Tumor", , drop = FALSE]

# ---- prepare plotting input ----
plot_df <- as.data.frame(t(cytosig_df))
plot_df$subtype <- clinical_tumor_df$LME_subtype

# ---- plot per cytokine ----
for (cyto_i in c("TGFB1", "GCSF", "TNFA", "IFNG")) {
  p <- ggplot(plot_df, aes(x = subtype, y = .data[[cyto_i]])) +
    geom_jitter(aes(fill = subtype), width = 0.2, alpha = 0.8,
                size = 0.75, stroke = 0.116, shape = 21) +
    geom_boxplot(outlier.shape = NA, size = 0.116, width = 0.3,
                 color = "black", fill = NA) +
    labs(x = "", y = "Signaling activity", title = cyto_i) +
    scale_fill_manual(values = subtype_col1) +
    theme(
      text = element_text(size = 5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line  = element_line(colour = "black", size = 0.116),
      axis.ticks = element_line(colour = "black", size = 0.116),
      axis.ticks.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.text.y  = element_text(colour = "black", size = 5),
      plot.title   = element_text(size = 5, hjust = 0.5),
      legend.position = "right",
      legend.title    = element_text(size = 5),
      legend.text     = element_text(size = 5)
    ) +
    ggpubr::stat_pwc(aes(group = subtype),
                     label.size = 5 * 25.4 / 72,
                     tip.length = 0,
                     size = 0.116,
                     hide.ns = TRUE,
                     vjust = 0.5,
                     method = "wilcox.test",
                     p.adjust.method = "BH",
                     label = "p.adj.signif") +
    coord_cartesian(clip = "off")
  
  ggsave(
    filename = file.path(fig_dir, paste0("fig1j.", cyto_i, ".pdf")),
    egg::set_panel_size(p, width = unit(1.5, "cm"), height = unit(1.5, "cm")),
    width = 8, height = 5, units = "cm"
  )
}