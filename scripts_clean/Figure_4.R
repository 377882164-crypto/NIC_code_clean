# ============================================================
# Figure 4 — Domain-Specific Motor Change (Pre–Post)
# From wide data → domain composites → change → plot
# ============================================================

library(tidyverse)

project_root <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw"
input_file <- file.path(project_root, "01_raw_data", "wide_dose_all3_cleaned_noHz.csv")
output_dir <- file.path(project_root, "06_figures", "Figure4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- read.csv(input_file)

data$group <- factor(
  data$group,
  levels = c(1, 2),
  labels = c("NeuroEPO", "Placebo")
)

fine_vars <- c(
  "fingertapr", "fingertapl",
  "handmover", "handmovel",
  "pronesupinerh", "pronesupinelh",
  "taprf", "taplf",
  "agilityrl", "agilityll"
)

rigidity_vars <- c(
  "rigidneck", "rigidupperr", "rigidupperl",
  "rigidlowerr", "rigidlowerl"
)

tremor_vars <- c(
  "posturaltremorrh", "posturaltremorlh",
  "actiontremorrh", "actiontremorlh",
  "restingtremorupperr", "restingtremorupperl",
  "restingtremorlowerr", "restingtremorlowerl",
  "restingtremorlip", "tremorpersistency"
)

agility_vars <- c(
  "chair", "walk", "freezerh",
  "posturalinstability", "postura", "movementspontaneity"
)

make_composite <- function(df, vars, name) {
  pre_cols <- paste0(vars, "_pre")
  post_cols <- paste0(vars, "_post")

  df %>%
    mutate(
      !!paste0(name, "_pre") := rowMeans(select(., all_of(pre_cols)), na.rm = TRUE),
      !!paste0(name, "_post") := rowMeans(select(., all_of(post_cols)), na.rm = TRUE)
    )
}

data <- data %>%
  make_composite(fine_vars, "fine") %>%
  make_composite(agility_vars, "agility") %>%
  make_composite(rigidity_vars, "rigidity") %>%
  make_composite(tremor_vars, "tremor")

data <- data %>%
  mutate(
    fine_change = fine_post - fine_pre,
    agility_change = agility_post - agility_pre,
    rigidity_change = rigidity_post - rigidity_pre,
    tremor_change = tremor_post - tremor_pre
  )

long_df <- data %>%
  select(
    group,
    fine_change,
    agility_change,
    rigidity_change,
    tremor_change
  ) %>%
  pivot_longer(
    cols = -group,
    names_to = "domain",
    values_to = "change"
  ) %>%
  mutate(
    domain = recode(
      domain,
      fine_change = "Fine motor",
      agility_change = "Agility",
      rigidity_change = "Rigidity",
      tremor_change = "Tremor"
    ),
    domain = factor(
      domain,
      levels = c("Fine motor", "Agility", "Rigidity", "Tremor")
    )
  )

summary_df <- long_df %>%
  group_by(group, domain) %>%
  summarise(
    n = sum(!is.na(change)),
    mean_change = mean(change, na.rm = TRUE),
    sd = sd(change, na.rm = TRUE),
    se = sd / sqrt(n),
    ci = 1.96 * se,
    ymin = mean_change - ci,
    ymax = mean_change + ci,
    .groups = "drop"
  )

COL_NEUROEPO <- "#F8766D"
COL_PLACEBO <- "#00BFC4"

p_fig4_final <- ggplot(
  summary_df,
  aes(
    x = domain,
    y = mean_change,
    color = group,
    group = group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey60",
    linewidth = 0.6
  ) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.2,
    linewidth = 0.9,
    position = position_dodge(width = 0.4)
  ) +
  scale_color_manual(
    values = c("NeuroEPO" = COL_NEUROEPO, "Placebo" = COL_PLACEBO)
  ) +
  labs(
    title = "Domain-Specific Motor Change (Pre-Post)",
    subtitle = "Mean change (post − pre) with 95% confidence intervals (±1.96 SE). Negative values indicate improvement.",
    x = NULL,
    y = "Change score (post − pre)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title.position = "panel",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave(
  file.path(output_dir, "Figure4_Domain_Level_Change.png"),
  p_fig4_final,
  width = 7.5,
  height = 4.8,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Figure4_Domain_Level_Change.pdf"),
  p_fig4_final,
  width = 7.5,
  height = 4.8,
  device = cairo_pdf,
  bg = "white"
)

message("Figure 4 exported (domain composites + pre-post change).")
