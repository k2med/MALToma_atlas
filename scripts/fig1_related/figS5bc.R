# ---- packages ----
library(survival)
library(survminer)
library(rms)
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

# ---- define variables ----
clinical_var_vec <- c(
  "Sex",
  "BIRC3_MALT_fusion",
  "Age_group",
  "Stage_group",
  "Site_group",
  "Subtype_group"
)

clinical_tumor_df[clinical_var_vec] <- lapply(
  clinical_tumor_df[clinical_var_vec],
  function(x) {
    if (is.character(x)) factor(x) else x
  }
)

# ---- univariate cox ----
unicox_results_df <- lapply(clinical_var_vec, function(var_name) {
  formula_obj <- as.formula(
    paste0("Surv(PFS_time, Progression_status) ~ ", var_name)
  )
  
  cox_fit <- coxph(formula_obj, data = clinical_tumor_df)
  cox_summary <- summary(cox_fit)
  
  data.frame(
    Variable   = var_name,
    Level      = rownames(cox_summary$coefficients),
    HR         = cox_summary$conf.int[, "exp(coef)"],
    Lower95CI  = cox_summary$conf.int[, "lower .95"],
    Upper95CI  = cox_summary$conf.int[, "upper .95"],
    P.value    = cox_summary$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

# ---- prepare univariate labels ----
unicox_results_df <- unicox_results_df %>%
  mutate(
    Label = case_when(
      Variable == "Site_group"    ~ gsub("Site_group", "Site: ", Level),
      Variable == "Age_group"     ~ gsub("Age_group", "Age: ", Level),
      Variable == "Subtype_group" ~ gsub("Subtype_group", "Subtype: ", Level),
      TRUE                        ~ Variable
    ),
    p_label = case_when(
      P.value < 0.001 ~ "***",
      P.value < 0.01  ~ "**",
      P.value < 0.05  ~ "*",
      TRUE            ~ NA_character_
    ),
    y = ifelse(HR > 1, Upper95CI + 0.5, Lower95CI - 0.5)
  )

unicox_order_vec <- c(
  "Site: Others",
  "Site: Stomach",
  "BIRC3_MALT_fusion",
  "Subtype: LME_2",
  "Stage_group",
  "Sex",
  "Age: >=60"
)

unicox_results_df <- unicox_results_df %>%
  mutate(Label = factor(Label, levels = unicox_order_vec))

# ---- forest plot function ----
save_forest_plot <- function(plot_df, x_var, out_name,
                             panel_width = 1.5,
                             panel_height = 2.0) {
  forest_plot <- ggplot(plot_df) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.232) +
    geom_linerange(
      aes(x = .data[[x_var]], ymin = Lower95CI, ymax = Upper95CI),
      linewidth = 0.232,
      show.legend = FALSE
    ) +
    geom_point(
      aes(x = .data[[x_var]], y = HR),
      size = 1,
      pch = 19,
      color = "#7f9db4"
    ) +
    geom_text(
      aes(x = .data[[x_var]], y = y, label = p_label),
      size = 5 / .pt,
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = c(0, 0.5)) +
    scale_y_log10(breaks = c(1, 10, 50), labels = c("1", "10", "50")) +
    labs(x = NULL, y = "Hazard ratio") +
    coord_flip() +
    theme_bw() +
    theme(
      text = element_text(size = 5),
      panel.grid.major = element_line(colour = "grey90", size = 0.116),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", size = 0.116),
      axis.line = element_line(colour = "black", size = 0.116),
      axis.ticks = element_line(colour = "black", size = 0.116),
      axis.text.x = element_text(colour = "black", size = 5),
      axis.text.y = element_text(colour = "black", size = 5),
      legend.position = "none",
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 5)
    )
  
  ggsave(
    filename = file.path(fig_dir, out_name),
    plot = egg::set_panel_size(
      forest_plot,
      width = grid::unit(panel_width, "cm"),
      height = grid::unit(panel_height, "cm")
    ),
    width = 6,
    height = 6,
    units = "cm"
  )
}

# ---- save univariate forest plot ----
save_forest_plot(
  plot_df = unicox_results_df,
  x_var = "Label",
  out_name = "figS5b_top.pdf"
)

# ---- multivariable cox ----
sig_var_vec <- unique(
  unicox_results_df$Variable[
    !is.na(unicox_results_df$P.value) & unicox_results_df$P.value < 0.05
  ]
)

if (length(sig_var_vec) == 0) {
  message("No variables with p < 0.05 in univariate Cox. Multivariable Cox skipped.")
} else {
  multicox_formula <- as.formula(
    paste0("Surv(PFS_time, Progression_status) ~ ", paste(sig_var_vec, collapse = " + "))
  )
  
  multicox_fit <- coxph(multicox_formula, data = clinical_tumor_df)
  multicox_summary <- summary(multicox_fit)
  
  multicox_results_df <- data.frame(
    Variable  = rownames(multicox_summary$conf.int),
    HR        = multicox_summary$conf.int[, "exp(coef)"],
    Lower95CI = multicox_summary$conf.int[, "lower .95"],
    Upper95CI = multicox_summary$conf.int[, "upper .95"],
    P.value   = multicox_summary$coefficients[, "Pr(>|z|)"],
    row.names = NULL,
    check.names = FALSE
  ) %>%
    mutate(
      p_label = case_when(
        P.value < 0.001 ~ "***",
        P.value < 0.01  ~ "**",
        P.value < 0.05  ~ "*",
        TRUE            ~ NA_character_
      ),
      y = ifelse(HR > 1, Upper95CI + 0.5, Lower95CI - 0.5)
    )
  
  write.table(
    multicox_results_df,
    file = "multicox.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  # ---- save multivariable forest plot ----
  save_forest_plot(
    plot_df = multicox_results_df,
    x_var = "Variable",
    out_name = "figS5b_bottom.pdf"
  )
  
  # ---- nomogram ----
  dd <- datadist(clinical_tumor_df)
  options(datadist = "dd")
  
  nomogram_fit <- cph(
    multicox_formula,
    x = TRUE,
    y = TRUE,
    surv = TRUE,
    data = clinical_tumor_df,
    time.inc = 60
  )
  
  surv_fun <- Survival(nomogram_fit)
  
  nomogram_obj <- nomogram(
    nomogram_fit,
    fun = list(
      function(x) surv_fun(12, x),
      function(x) surv_fun(36, x),
      function(x) surv_fun(60, x)
    ),
    lp = FALSE,
    funlabel = c("1-year PFS", "3-year PFS", "5-year PFS"),
    maxscale = 100,
    fun.at = c(0.99, 0.9, 0.8, 0.6, 0.4, 0.2, 0.1)
  )
  
  pdf(file = "figS5c.pdf", height = 5, width = 10)
  plot(nomogram_obj)
  dev.off()
}
