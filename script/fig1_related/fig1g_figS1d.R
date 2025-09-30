# ---- packages ----
library(ggplot2)
library(ggpubr)
library(dplyr)
library(egg)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/bulk_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

bulk_clinical_path <- file.path(data_dir, "bulk_clinical_metadata.txt")
source(file.path(src_dir, "custom_colors.R"))

# ---- load data ----
clinical_df <- read.table(
  bulk_clinical_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE,
  na.strings = c("", "NA")
)
clinical_tumor_df <- clinical_df[clinical_df$Group == "Tumor", , drop = FALSE]

# ---- derive grouped variables (keep original column name casing) ----
clinical_tumor_df <- clinical_tumor_df %>%
  mutate(
    Age_group = case_when(Age < 60 ~ "<60", Age >= 60 ~ ">=60"),
    Stage_group = case_when(
      Ann_Arbor_stage %in% c("I", "II", "III") ~ "Stage I_III",
      Ann_Arbor_stage %in% c("IV") ~ "Stage IV",
      TRUE ~ NA_character_
    ),
    Site_group = if_else(Site %in% c("Ocular_adnexa", "Stomach"), Site, "Others")
  )

# ---- helper: auto test (chisq vs fisher) and stacked proportional barplot ----
plot_prop_bar <- function(df, var, lme_col = "LME_subtype",
                          fill_values, out_name,
                          x_text_angle = 60,
                          x_levels = NULL) {
  
  plot_df <- df[!is.na(df[[var]]) & !is.na(df[[lme_col]]), , drop = FALSE]
  
  if (!is.null(x_levels) && var %in% names(plot_df)) {
    plot_df[[var]] <- factor(plot_df[[var]], levels = x_levels)
  }
  
  tab <- table(plot_df[[var]], plot_df[[lme_col]])
  if (any(tab < 5)) {
    test_result <- fisher.test(tab)
  } else {
    test_result <- suppressWarnings(chisq.test(tab))
  }
  pval <- signif(test_result$p.value, 3)
  
  p <- ggplot(plot_df, aes(x = .data[[lme_col]], fill = .data[[var]])) +
    geom_bar(position = "fill", width = 0.6) +
    labs(
      title = paste0(var, " distribution by ", lme_col),
      subtitle = paste0("P = ", pval),
      x = "",
      y = "Proportion (%)",
      fill = var
    ) +
    scale_fill_manual(values = fill_values) +
    theme_minimal() +
    theme(
      text = element_text(size = 5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", size = 0.116),
      axis.ticks = element_line(colour = "black", size = 0.116),
      axis.text.x = element_text(colour = "black", size = 5,
                                 angle = x_text_angle, vjust = 1, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 5),
      plot.title = element_text(size = 5, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 5)
    ) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.5),
                       labels = function(x) x * 100)
  
  ggsave(
    filename = file.path(fig_dir, out_name),
    plot = egg::set_panel_size(p, width = unit(1.5, "cm"), height = unit(2, "cm")),
    width = 8, height = 5, units = "cm"
  )
}

# ---- plots ----
# 1) Sex
plot_prop_bar(
  df = clinical_tumor_df,
  var = "Sex",
  fill_values = sex_col,
  out_name = "fig1g_right.pdf"
)

# 2) Site_group (ordered)
plot_prop_bar(
  df = clinical_tumor_df,
  var = "Site_group",
  fill_values = site_col,
  out_name = "fig1g_left.pdf",
  x_levels = c("Ocular_adnexa", "Stomach", "Others")
)

# 3) Stage_group (ordered)
plot_prop_bar(
  df = clinical_tumor_df,
  var = "Stage_group",
  fill_values = stage_group_col,
  out_name = "figS1d_left.pdf",
  x_levels = c("Stage I_III", "Stage IV")
)

# 4) Age_group (ordered)
plot_prop_bar(
  df = clinical_tumor_df,
  var = "Age_group",
  fill_values = age_group_col,
  out_name = "figS1d_right.pdf",
  x_levels = c("<60", ">=60")
)