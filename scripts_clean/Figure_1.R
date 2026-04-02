# ============================================================
# Make Figure 1B heatmap from GLMM interaction CSV
# Output: PNG (600 dpi) + PDF
# FIXED: subtitle wrapping to avoid truncation
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(RColorBrewer)
  library(stringr)
})

# ---------------------------
# Paths (EDIT THESE)
# ---------------------------
path_glmm_topz <- "D:/r/nic/stats_GLMM_topz_interaction.csv"
out_dir        <- "D:/r/nic/Figure1_remake"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Export parameters
# ---------------------------
W_IN <- 8.6
H_IN <- 5.0
DPI  <- 600

# ---------------------------
# Orders
# ---------------------------
bands_order   <- c("Delta", "Theta", "LowAlpha", "HighAlpha", "Beta")
regions_order <- c("Frontal", "Central", "Parietal", "Occipital", "Temporal")

# ---------------------------
# Helpers
# ---------------------------
save_png_pdf <- function(p, filename_no_ext, width, height, dpi = 600) {
  ggsave(
    filename = file.path(out_dir, paste0(filename_no_ext, ".png")),
    plot = p, width = width, height = height, units = "in", dpi = dpi,
    bg = "white", limitsize = FALSE
  )
  ggsave(
    filename = file.path(out_dir, paste0(filename_no_ext, ".pdf")),
    plot = p, width = width, height = height, units = "in",
    device = cairo_pdf, bg = "white", limitsize = FALSE
  )
}

star_from_padj <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

legend_breaks <- function(lim_z) {
  br <- c(-2, -1, 0, 1, 2)
  br[br >= -lim_z & br <= lim_z]
}

prep_interaction_df <- function(path_csv) {
  df <- readr::read_csv(path_csv, show_col_types = FALSE) %>%
    filter(term == "GroupNeuroEPO:Timepost") %>%
    mutate(
      Band      = factor(as.character(Band), levels = bands_order),
      Region    = factor(as.character(Region), levels = regions_order),
      estimate  = suppressWarnings(as.numeric(estimate)),
      std.error = if ("std.error" %in% names(.)) suppressWarnings(as.numeric(.data[["std.error"]])) else NA_real_,
      statistic = if ("statistic" %in% names(.)) suppressWarnings(as.numeric(statistic)) else NA_real_,
      p_adj     = suppressWarnings(as.numeric(p_adj)),
      z = dplyr::case_when(
        is.finite(statistic) ~ statistic,
        is.finite(std.error) & std.error > 0 ~ estimate / std.error,
        TRUE ~ estimate
      ),
      star = if_else(is.finite(p_adj) & p_adj < 0.05, star_from_padj(p_adj), "")
    )
  
  df2 <- df %>%
    filter(Band %in% bands_order, Region %in% regions_order) %>%
    arrange(Band, Region)
  
  # Must be exactly 25 cells
  expected <- length(bands_order) * length(regions_order)
  if (nrow(df2) != expected) {
    stop(
      "Expected ", expected, " rows (5 bands × 5 regions) but got ", nrow(df2),
      " in file: ", path_csv
    )
  }
  
  # Must be unique per Band×Region
  dup <- df2 %>% count(Band, Region) %>% filter(n != 1)
  if (nrow(dup) > 0) {
    stop("Band×Region not unique in: ", path_csv)
  }
  
  df2
}

plot_interaction_heatmap <- function(df, title, subtitle, lim_z) {
  pal  <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  brks <- legend_breaks(lim_z)
  
  subtitle_wrapped <- stringr::str_wrap(subtitle, width = 95)
  
  ggplot(df, aes(x = Region, y = Band, fill = z)) +
    geom_tile(color = "white", linewidth = 0.6, alpha = 0.7, na.rm = FALSE) +
    geom_text(aes(label = star), size = 5, color = "black", na.rm = TRUE) +
    scale_fill_gradientn(
      colours  = pal,
      values   = scales::rescale(c(-lim_z, 0, lim_z)),
      limits   = c(-lim_z, lim_z),
      breaks   = brks,
      labels   = brks,
      oob      = scales::squish,
      na.value = "grey90",
      name     = NULL
    ) +
    labs(
      title = title,
      subtitle = subtitle_wrapped,
      x = "Region",
      y = "Band"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 11) +
    theme(
      plot.title.position = "panel",
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(
        size = 7.8,
        lineheight = 1.05,
        hjust = 0,
        margin = margin(b = 6)
      ),
      axis.title = element_text(size = 11),
      axis.text  = element_text(size = 10, colour = "black"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "right",
      plot.margin = margin(12, 18, 32, 12)
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_discrete(expand = expansion(mult = c(0.02, 0.02)))
}

# ---------------------------
# Load + prep
# ---------------------------
dfB <- prep_interaction_df(path_glmm_topz)

# Color limits (95th percentile of |z|)
q_level <- 0.95
lim_z <- as.numeric(quantile(abs(dfB$z), probs = q_level, na.rm = TRUE))
if (!is.finite(lim_z) || lim_z <= 0) lim_z <- max(abs(dfB$z), na.rm = TRUE)
if (!is.finite(lim_z) || lim_z <= 0) lim_z <- 1

subtitle_text <- paste0(
  "Binomial GLMM (GroupNeuroEPO × Timepost). ",
  "Fill: Wald test statistic (symmetric limits at ", q_level * 100,
  "th percentile of |statistic|). ",
  "Stars: BH-FDR p_adj < 0.05."
)

# ---------------------------
# Make plot
# ---------------------------
pB <- plot_interaction_heatmap(
  dfB,
  title = "Top-Z hub rate: Group × Time interaction",
  subtitle = subtitle_text,
  lim_z = lim_z
)

# ---------------------------
# Save
# ---------------------------
save_png_pdf(
  pB,
  "Figure1B_topz_interaction_heatmap",
  width = W_IN,
  height = H_IN,
  dpi = DPI
)

message("Done. Files saved to: ", normalizePath(out_dir))