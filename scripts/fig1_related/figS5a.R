# ---- packages ----
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(ggpubr)

# ---- paths & sources ----
data_dir  <- "../../data/bulk_rnaseq"
src_dir   <- "../source"
fig_dir   <- "."

bulk_clinical_path <- file.path(data_dir, "bulk_clinical_metadata.txt")
source(file.path(src_dir, "custom_colors.R")) 

# ---- data: clinical (tumor only) ----
clinical_df <- read.table(
  bulk_clinical_path, header = TRUE, sep = "\t", row.names = 1, na.strings = c("", "NA")
)

clinical_tumor_df <- clinical_df[clinical_df$Group == "Tumor", , drop = FALSE] %>%
  mutate(
    Age_group = case_when(Age < 60 ~ "<60", Age >= 60 ~ ">=60"),
    Stage_group = case_when(
      Ann_Arbor_stage %in% c("I", "II", "III") ~ "Stage I_III",
      Ann_Arbor_stage %in% c("IV") ~ "Stage IV",
      TRUE ~ NA_character_
    ),
    Site_group = if_else(Site %in% c("Ocular_adnexa", "Stomach"), Site, "Others"),
    Subtype_group = if_else(LME_subtype %in% c("LME_1", "LME_3"), "LME_1/3", LME_subtype)
  )

# ---- helper: survival plot by group ----
plot_survival_by_group <- function(clinical_data, group_var, palette_vector,
                                   title_text = "PFS",
                                   show_risktable = TRUE,
                                   show_cumevents = TRUE,
                                   risk_table_height = 0.3,
                                   time_col = "PFS_time", status_col = "Progression_status",
                                   outfile_prefix = NA) {
  
  keep_cols <- c(status_col, time_col, group_var)
  stopifnot(all(keep_cols %in% colnames(clinical_data)))
  
  df <- clinical_data[, keep_cols, drop = FALSE]
  colnames(df) <- c("status", "time", "group")
  
  # drop NAs
  df <- df %>% filter(!is.na(status), !is.na(time), !is.na(group))
  
  surv_fit <- survfit(Surv(time, status) ~ group, data = df)
  
  surv_plot <- ggsurvplot(
    fit = surv_fit,
    data = df,
    conf.int = FALSE,
    pval = TRUE,
    surv.median.line = "hv",
    title = title_text,
    xlab = "Time (month)",
    palette = unname(palette_vector[levels(factor(df$group))]),
    legend = "top",
    legend.title = group_var,
    legend.labs = levels(df$group),
    risk.table = show_risktable,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    cumevents = show_cumevents,
    cumevents.y.text.col = TRUE,
    cumevents.y.text = FALSE,
  )
  
  final_plot <- ggpubr::ggarrange(
    surv_plot$plot,
    surv_plot$table,
    surv_plot$cumevents,
    ncol = 1,
    align = "v",
    heights = c(0.7, risk_table_height, risk_table_height)
  )
  
  pdf(file.path(fig_dir, paste0(outfile_prefix, ".", group_var, ".pdf")),
      width = 4, height = 6, onefile = FALSE)
  print(final_plot)
  dev.off()
}

# ---- run plots ----
plot_survival_by_group(clinical_tumor_df, "Stage_group", stage_group_col, outfile_prefix = "figS5a",
                       show_risktable = TRUE, risk_table_height = 0.3)
plot_survival_by_group(clinical_tumor_df, "BIRC3_MALT_fusion", fusion_col, outfile_prefix = "figS5a",
                       show_risktable = TRUE, risk_table_height = 0.3)
plot_survival_by_group(clinical_tumor_df, "TNFAIP3_mutation", mutation_col, outfile_prefix = "figS5a",
                       show_risktable = TRUE, risk_table_height = 0.3)
plot_survival_by_group(clinical_tumor_df, "Site_group", site_col, outfile_prefix = "figS5a",
                       show_risktable = TRUE, risk_table_height = 0.35)
plot_survival_by_group(clinical_tumor_df, "Age_group", age_group_col, outfile_prefix = "figS5a",
                       show_risktable = TRUE, risk_table_height = 0.3)
plot_survival_by_group(clinical_tumor_df, "Sex", sex_col, outfile_prefix = "figS5a",
                       show_risktable = TRUE, risk_table_height = 0.3)
plot_survival_by_group(clinical_tumor_df, "LME_subtype", subtype_col2, outfile_prefix = "figS5a",
                       show_risktable = TRUE, risk_table_height = 0.35)
plot_survival_by_group(clinical_tumor_df, "Subtype_group", subtype_col3, outfile_prefix = "fig1i",
                       show_risktable = TRUE, risk_table_height = 0.3)
