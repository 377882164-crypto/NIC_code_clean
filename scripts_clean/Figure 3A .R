# =============================
# Figure 3A — Forest plot (CORRECTED)
# Top-Z rate; GENUINE "both" criterion; standardized indirect effects
# Improvement direction using variable dictionary
# =============================

library(tidyverse)
library(readxl)
library(janitor)
library(RColorBrewer)
library(scales)
library(stringr)
library(grid)

# ---------- PATHS (EDIT THESE) ----------
path_topz_both <- "D:/r/nic/mediation_results_main_split/main_split__TOPZ_BOTH_Genuine_key_results.csv"
path_dict      <- "D:/r/nic修改/csv/NeuroEPO PD Motor & Cognition Variable Dictionary_notitle.xlsx"

out_png <- "D:/r/nic/3a/Figure3A_Mediation_Forest_Dose_STD_TopZ_BOTH_Genuine.png"

# ---------- HELPERS ----------
canon_abbrev <- function(x) {
  x %>%
    sub("_(post|pre)$", "", ., ignore.case = TRUE) %>%
    sub("_a\\d+_\\d+$", "", ., ignore.case = TRUE) %>%
    sub("_a\\d+$", "", ., ignore.case = TRUE)
}

# Shortened labels for publication
short_names <- c(
  "pronesupinelh"       = "Pronation–supination (Left)",
  "rigidlowerl"         = "Rigidity (Left lower limb)",
  "rigidupperl"         = "Rigidity (Left upper limb)",
  "taplf"               = "Toe tapping (Left)",
  "restingtremorlowerr" = "Resting tremor (Right lower limb)",
  "speech"              = "Speech",
  "wmSequence"          = "WM sequencing (accuracy)",
  "tiempototal2"        = "Stroop word-reading time"
)

# ---------- CREATE OUTPUT DIR ----------
out_dir <- dirname(out_png)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ---------- READ DATA ----------
dat <- readr::read_csv(path_topz_both, show_col_types = FALSE) %>%
  filter(mediator_metric == "topz_rate") %>%
  # Safety: keep only rows with finite standardized CI
  filter(is.finite(ind_est_std), is.finite(ind_lwr_std), is.finite(ind_upr_std))

if (nrow(dat) == 0) {
  stop("❌ No results found in input file. Check TOPZ_BOTH_Genuine_key_results.csv generation.")
}

cat("✓ Loaded", nrow(dat), "TopZ mediation results meeting BOTH criterion\n")

dict <- readxl::read_excel(path_dict, sheet = 1) %>%
  janitor::clean_names() %>%
  transmute(
    abbrev    = str_trim(abbreviation),
    direction = str_trim(direction),
    full_name = str_trim(full_english_name),
    catalog   = str_trim(catalog)
  )

# ---------- DICTIONARY JOIN ----------
plot_data <- dat %>%
  mutate(y_abbrev = canon_abbrev(y_post)) %>%
  left_join(dict, by = c("y_abbrev" = "abbrev"))

# Fail fast if dictionary missing
if (any(is.na(plot_data$direction))) {
  missing <- plot_data %>% filter(is.na(direction)) %>% distinct(y_post, y_abbrev)
  stop(
    "Missing direction in dictionary for:\n",
    paste0(missing$y_post, " (canon: ", missing$y_abbrev, ")", collapse = "\n")
  )
}

# ---------- IMPROVEMENT DIRECTION FLIP ----------
plot_data <- plot_data %>%
  mutate(
    lower_better = str_to_lower(direction) %>% str_detect("^lower is better"),
    ind_est_adj  = if_else(lower_better, -ind_est_std, ind_est_std),
    ind_lwr_adj  = if_else(lower_better, -ind_upr_std, ind_lwr_std),
    ind_upr_adj  = if_else(lower_better, -ind_lwr_std, ind_upr_std)
  )

# ---------- LABELS ----------
plot_data <- plot_data %>%
  mutate(
    outcome_label = coalesce(
      unname(short_names[y_abbrev]),
      na_if(full_name, ""),
      y_abbrev
    ),
    band_region = paste0(mediator_band, " ", mediator_region, " (Top‑Z)"),
    row_label   = paste0(mediator_band, " ", mediator_region, " → ", outcome_label)
  ) %>%
  arrange(mediator_band, mediator_region, outcome_label)

# Order so first row appears at top
plot_data$row_label <- factor(plot_data$row_label, levels = rev(plot_data$row_label))

# ---------- COLORS ----------
band_levels <- sort(unique(plot_data$band_region))

pal <- setNames(
  brewer.pal(max(3, min(length(band_levels), 8)), "Dark2")[seq_along(band_levels)],
  band_levels
)

# Symmetric x-limits
xmax <- max(abs(c(plot_data$ind_lwr_adj, plot_data$ind_upr_adj)), na.rm = TRUE) * 1.1
xlim <- c(-xmax, xmax)

# Wrapped title to avoid truncation
plot_title <- stringr::str_wrap(
  "Dose-Dependent Mediation of Motor Change via Top‑Z Rate Reorganization",
  width = 55
)

# ---------- PLOT ----------
p_forest <- ggplot(plot_data, aes(x = ind_est_adj, y = row_label, color = band_region)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.9) +
  geom_errorbarh(
    aes(xmin = ind_lwr_adj, xmax = ind_upr_adj),
    height = 0.25,
    linewidth = 1.1
  ) +
  geom_point(size = 3.6, shape = 19) +
  scale_color_manual(values = pal, drop = FALSE) +
  scale_x_continuous(
    limits = xlim,
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Standardized Indirect Effect (Improvement Direction)",
    y = NULL,
    title = plot_title
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 15,
      lineheight = 1.1,
      margin = margin(b = 10)
    ),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    plot.margin = margin(t = 15, r = 30, b = 15, l = 15)
  )

ggsave(
  filename = out_png,
  plot = p_forest,
  width = 12.5,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

cat("✓ Figure saved:", out_png, "\n")
cat("✓ Total pathways plotted:", nrow(plot_data), "\n")

p_forest