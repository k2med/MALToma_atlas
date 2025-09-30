# ---- packages ----
library(GSVA)
library(GSEABase)
library(ggplot2)
library(survival)
library(survminer)
library(dplyr)
library(grid)

# ---- paths & sources ----
data_dir <- "../../data/bulk_rnaseq"
src_dir  <- "../source"
fig_dir  <- "."

bulk_expr_path     <- file.path(data_dir, "bulk_expression_matrix.txt")
bulk_clinical_path <- file.path(data_dir, "bulk_clinical_metadata.txt")
gmt_path           <- file.path(src_dir, "In_house_LME_signatures.gmt")

source(file.path(src_dir, "custom_colors.R"))

# ---- clinical data (tumor only) ----
clinical_df <- read.table(
  bulk_clinical_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)
clinical_tumor_df <- clinical_df[clinical_df$Group == "Tumor", , drop = FALSE]

# ---- load expression & GSVA ----
bulk_expr_mat <- read.table(
  bulk_expr_path, header = TRUE, row.names = 1, check.names = FALSE
)
# align expression to tumor samples only
bulk_expr_mat <- bulk_expr_mat[, rownames(clinical_tumor_df), drop = FALSE]
bulk_expr_mat <- as.matrix(log2(bulk_expr_mat + 1))

gene_sets <- getGmt(gmt_path, sep = "\t", geneIdType = SymbolIdentifier())
gsva_score_mat <- gsva(bulk_expr_mat, gene_sets, method = "gsva", kcdf = "Gaussian")

# ---- B_Mal signature split & PFS KM ----
# append B_Mal signature score to clinical data
clinical_tumor_df <- cbind(
  clinical_tumor_df,
  t(gsva_score_mat["B_Mal", rownames(clinical_tumor_df), drop = FALSE])
)

clinical_tumor_df <- clinical_tumor_df %>%
  mutate(B_Mal_group = if_else(B_Mal > median(B_Mal, na.rm = TRUE), "High", "Low"))

pfs_fit <- survfit(Surv(PFS_time, Progression_status) ~ B_Mal_group, data = clinical_tumor_df)

pfs_plot <- ggsurvplot(
  pfs_fit,
  data = clinical_tumor_df,
  conf.int = FALSE,
  pval = TRUE,
  surv.median.line = "hv",
  xlab = "Time (month)",
  palette = c("#d5201b", "#3777ac"),
  legend = "top",
  legend.title = "B_Mal signature",
  legend.labs = c("High", "Low"),
  risk.table = FALSE
)

pdf(file = file.path(fig_dir, "fig3e.pdf"), width = 4, height = 4, onefile = F)
print(pfs_plot)
dev.off()

# ---- re-run GSVA on MP gene sets from CSV ----
gene_sets_df <- read.csv(
  file.path(src_dir, "MP_list.csv"),
  header = TRUE, check.names = FALSE, stringsAsFactors = FALSE
)
gene_sets_list <- as.list(gene_sets_df)

gsva_score_mat <- gsva(bulk_expr_mat, gene_sets_list, method = "gsva", kcdf = "Gaussian")

# append all MP scores to clinical data (samples as rows)
clinical_tumor_df <- cbind(
  clinical_tumor_df,
  t(gsva_score_mat[, rownames(clinical_tumor_df), drop = FALSE])
)

# ---- univariate Cox for each MP ----
cox_results <- lapply(names(gene_sets_list), function(mp_name) {
  clinical_input <- clinical_tumor_df[, c("Progression_status", "PFS_time", mp_name)]
  colnames(clinical_input)[ncol(clinical_input)] <- "mp_score"
  
  cox_model <- coxph(Surv(PFS_time, Progression_status) ~ mp_score, data = clinical_input)
  sm <- summary(cox_model)
  
  data.frame(
    Variable   = mp_name,
    HR         = unname(sm$coefficients[1, "exp(coef)"]),
    Lower95CI  = unname(sm$conf.int[1, "lower .95"]),
    Upper95CI  = unname(sm$conf.int[1, "upper .95"]),
    P.value    = unname(sm$coefficients[1, "Pr(>|z|)"]),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

# adjust p-values and prepare labels/positions
cox_results <- cox_results %>%
  mutate(
    p.adjusted = p.adjust(P.value, method = "BH"),
    p.adjusted_label = case_when(
      p.adjusted < 0.001 ~ "***",
      p.adjusted < 0.01  ~ "**",
      p.adjusted < 0.05  ~ "*",
      TRUE               ~ NA_character_
    ),
    y = ifelse(HR > 1, Upper95CI + 0.5, Lower95CI - 0.5)
  )

# ---- forest plot ----
cox_results$Variable <- factor(
  cox_results$Variable,
  levels = cox_results$Variable[order(cox_results$HR, decreasing = FALSE)]
)

forest_plot <- ggplot(cox_results) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.232) +
  geom_linerange(aes(x = Variable, ymin = Lower95CI, ymax = Upper95CI),
                 linewidth = 0.232, show.legend = FALSE) +
  geom_point(aes(x = Variable, y = HR), size = 1, pch = 19, color = "#7f9db4") +
  geom_text(aes(x = Variable, y = y, label = p.adjusted_label),
            size = 5 / .pt, show.legend = FALSE) +
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
    legend.text  = element_text(size = 5)
  )

ggsave(
  filename = file.path(fig_dir, "fig3f.pdf"),
  plot = egg::set_panel_size(
    forest_plot,
    width  = grid::unit(1.5, "cm"),
    height = grid::unit(2.0, "cm")
  ),
  width = 6, height = 6, units = "cm"
)