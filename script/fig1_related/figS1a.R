# ---- packages ----
library(dplyr)
library(ggplot2)
library(gg.gap)
library(egg)
library(grid)

# ---- paths ----
data_dir  <- "../../data/bulk_rnaseq"
src_dir   <- "../source"
fig_dir   <- "."

bulk_clinical_path <- file.path(data_dir, "bulk_clinical_metadata.txt")

source(file.path(src_dir, "custom_colors.R"))

# ---- load data ----
clinical_df <- read.table(
  bulk_clinical_path, header = TRUE, sep = "\t",
  row.names = 1, check.names = FALSE
)

# ---- prepare counts ----
tumor_counts_df <- clinical_df %>%
  filter(Group == "Tumor") %>%
  count(Site) %>%
  mutate(Group = "Tumor")

normal_count_df <- clinical_df %>%
  filter(Group == "Normal") %>%
  summarise(n = n()) %>%
  mutate(Site = "Normal", Group = "Normal")

plot_df <- bind_rows(tumor_counts_df, normal_count_df)

# ---- site order ----
site_order_vec <- tumor_counts_df %>%
  arrange(desc(n)) %>%
  pull(Site) %>%
  c("Normal")

plot_df$Site <- factor(plot_df$Site, levels = site_order_vec)

# ---- barplot ----
barplot_base <- ggplot(plot_df, aes(x = Site, y = n, fill = Group)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = group_col) +
  labs(x = "", y = "Count") +
  theme_classic() +
  theme(
    text             = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line        = element_line(colour = "black", size = 0.116),
    axis.ticks       = element_line(colour = "black", size = 0.116),
    axis.text.x      = element_text(colour = "black", size = 10,
                                    angle = 60, vjust = 1, hjust = 1),
    axis.text.y      = element_text(colour = "black", size = 10),
    plot.title       = element_text(size = 10, hjust = 0.5),
    legend.position  = "right",
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 10)
  )

barplot_gap <- gg.gap(
  plot        = barplot_base,
  segments    = c(10, 20),
  tick_width  = c(5, 10),
  rel_heights = c(1, 0, 1),
  ylim        = c(0, max(plot_df$n))
)

# ---- save figure ----
ggsave(
  filename = file.path(fig_dir, "figS1a.pdf"),
  egg::set_panel_size(barplot_gap,
                      width  = unit(10, "cm"),
                      height = unit(10, "cm")),
  width = 16, height = 16, units = "cm"
)
