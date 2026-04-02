################################################################################
# Figure 2 (clean mainline)
# A) stats_GLMM_topz_interaction.csv (LowAlpha)
# B) electrode topography from signal_activity_nodelevel_all.csv
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(patchwork)
})

# =========================
# Paths
# =========================
project_root <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw"
stats_file <- file.path(project_root, "05_results", "clean_topz_integrated", "topz_results", "stats_GLMM_topz_interaction.csv")
nodelevel_file <- file.path(project_root, "05_results", "clean_topz_integrated", "signal_activity_nodelevel_all.csv")
locs_file <- file.path(project_root, "02_auxiliary", "Standard-10-20-Cap18(1).locs")

output_dir <- file.path(project_root, "06_figures", "Figure2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

TOP_PROP <- 0.10

# =========================
# Helpers
# =========================
stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

std_group <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    tolower(x) %in% c("placebo", "pla") ~ "Placebo",
    tolower(x) %in% c("neuroepo", "neuro_epo", "neuro") ~ "NeuroEPO",
    TRUE ~ x
  )
}

std_time <- function(x) {
  x <- trimws(as.character(x))
  x2 <- tolower(x)
  dplyr::case_when(
    x2 %in% c("pre", "baseline") ~ "pre",
    x2 %in% c("post", "followup", "follow-up") ~ "post",
    TRUE ~ x2
  )
}

std_ch <- function(x) {
  x <- trimws(as.character(x))
  x <- toupper(x)
  x <- gsub("[^A-Z0-9]", "", x)
  x
}

is_lowalpha <- function(band_str) {
  x <- tolower(trimws(as.character(band_str)))
  grepl("low\\s*-?\\s*alpha", x)
}

# =========================
# Panel A
# =========================
create_panel_a <- function(stats_file) {
  stats_data <- readr::read_csv(stats_file, show_col_types = FALSE)

  need <- c("Band", "Region", "term", "estimate", "std.error", "p_adj")
  miss <- setdiff(need, names(stats_data))
  if (length(miss) > 0) stop("stats_file missing columns: ", paste(miss, collapse = ", "))

  lowa <- stats_data %>%
    filter(Band == "LowAlpha", term == "GroupNeuroEPO:Timepost") %>%
    mutate(
      Region = factor(as.character(Region),
                      levels = c("Frontal", "Central", "Parietal", "Occipital", "Temporal")),
      sig_label = stars(p_adj),
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error
    ) %>%
    arrange(Region)

  if (nrow(lowa) == 0) stop("No rows for LowAlpha interaction in stats_file.")

  y_lim <- range(c(lowa$ci_lower, lowa$ci_upper), na.rm = TRUE)
  pad <- 0.10 * diff(y_lim)
  if (!is.finite(pad) || pad == 0) pad <- 0.5
  y_lim <- c(y_lim[1] - pad, y_lim[2] + pad)

  ggplot(lowa, aes(x = Region, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray45", linewidth = 0.6) +
    geom_col(width = 0.72, fill = "#0D47A1", color = "black", linewidth = 0.35, alpha = 0.85) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.22, linewidth = 0.55) +
    geom_text(
      aes(label = sig_label, y = ifelse(estimate >= 0, ci_upper, ci_lower)),
      vjust = ifelse(lowa$estimate >= 0, -0.7, 1.6),
      size = 6, fontface = "bold", color = "#B71C1C"
    ) +
    coord_cartesian(ylim = y_lim) +
    labs(
      title = "A. Top-Z hub rate interaction in LowAlpha",
      x = "Brain region",
      y = "GLMM interaction estimate (log-odds)\n(GroupNeuroEPO × Timepost)"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# =========================
# Panel B
# =========================
create_panel_b <- function(nodelevel_file, locs_file, top_prop = 0.10,
                           alpha_fdr = 0.05,
                           test_method = c("welch", "wilcox")) {

  test_method <- match.arg(test_method)

  node_all <- readr::read_csv(nodelevel_file, show_col_types = FALSE)

  need <- c("Group", "ID", "Time", "Channel", "Band", "Value_z")
  miss <- setdiff(need, names(node_all))
  if (length(miss) > 0) stop("nodelevel_file missing columns: ", paste(miss, collapse = ", "))

  lowa_nodes <- node_all %>%
    transmute(
      Group = std_group(Group),
      ID = as.character(ID),
      Time = std_time(Time),
      Channel = std_ch(Channel),
      Band = trimws(as.character(Band)),
      Value_z = as.numeric(Value_z)
    ) %>%
    filter(
      Group %in% c("Placebo", "NeuroEPO"),
      Time %in% c("pre", "post"),
      is_lowalpha(Band),
      is.finite(Value_z),
      !is.na(Channel), Channel != ""
    )

  if (nrow(lowa_nodes) == 0) stop("No LowAlpha rows after filtering (check Band/Time/Group coding).")

  lowa_flag <- lowa_nodes %>%
    group_by(Group, ID, Time) %>%
    mutate(
      n_nodes = n(),
      k = pmax(1L, ceiling(n_nodes * top_prop)),
      rk = rank(abs(Value_z), ties.method = "first"),
      is_topz = as.integer(rk > (n_nodes - k))
    ) %>%
    ungroup()

  ch_rate <- lowa_flag %>%
    group_by(Group, ID, Time, Channel) %>%
    summarise(topz_rate_ch = mean(is_topz, na.rm = TRUE), .groups = "drop")

  ch_wide <- ch_rate %>% pivot_wider(names_from = Time, values_from = topz_rate_ch)
  if (!("pre" %in% names(ch_wide)) || !("post" %in% names(ch_wide))) {
    stop("After pivot_wider, missing pre or post column (check Time coding).")
  }

  ch_delta <- ch_wide %>%
    filter(is.finite(pre), is.finite(post)) %>%
    mutate(delta = post - pre)

  if (nrow(ch_delta) == 0) stop("No paired pre/post data; cannot compute change.")

  gmean <- ch_delta %>%
    group_by(Group, Channel) %>%
    summarise(
      mean_delta = mean(delta, na.rm = TRUE),
      n_ID = n_distinct(ID),
      .groups = "drop"
    )

  g_ne <- gmean %>%
    filter(Group == "NeuroEPO") %>%
    select(Channel, mean_delta_ne = mean_delta, n_ID_ne = n_ID)

  g_pl <- gmean %>%
    filter(Group == "Placebo") %>%
    select(Channel, mean_delta_pl = mean_delta, n_ID_pl = n_ID)

  ch_diff <- g_ne %>%
    inner_join(g_pl, by = "Channel") %>%
    mutate(
      diff_change = mean_delta_ne - mean_delta_pl,
      abs_change = abs(diff_change)
    )

  if (nrow(ch_diff) == 0) stop("No channels have paired data in BOTH groups.")

  ch_p <- ch_delta %>%
    group_by(Channel) %>%
    summarise(
      p = {
        d <- dplyr::cur_data()
        if (length(unique(d$Group)) < 2) {
          NA_real_
        } else if (test_method == "welch") {
          tryCatch(stats::t.test(delta ~ Group, data = d)$p.value,
                   error = function(e) NA_real_)
        } else {
          tryCatch(stats::wilcox.test(delta ~ Group, data = d, exact = FALSE)$p.value,
                   error = function(e) NA_real_)
        }
      },
      n_ne = sum(Group == "NeuroEPO"),
      n_pl = sum(Group == "Placebo"),
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      sig_fdr = is.finite(p_adj) & (p_adj < alpha_fdr)
    )

  ch_diff2 <- ch_diff %>%
    left_join(ch_p, by = "Channel")

  locs_all <- read.table(locs_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(locs_all) <- c("num", "angle_deg", "radius", "label")

  locs_all <- locs_all %>%
    mutate(
      label_std = std_ch(label),
      theta = (90 - angle_deg) * pi / 180,
      x = radius * cos(theta),
      y = radius * sin(theta)
    )

  cz_idx <- which(locs_all$label_std == "CZ")
  if (length(cz_idx) > 0) {
    locs_all$x[cz_idx] <- 0
    locs_all$y[cz_idx] <- 0
  }
  radius_max <- max(locs_all$radius, na.rm = TRUE)
  if (!is.finite(radius_max) || radius_max <= 0) stop("Invalid radius_max from locs_file.")

  topo <- locs_all %>%
    inner_join(ch_diff2, by = c("label_std" = "Channel"))

  if (nrow(topo) == 0) stop("No matching channels between loc labels and computed channel diffs.")

  theta_circle <- seq(0, 2 * pi, length.out = 400)
  circle <- data.frame(x = radius_max * cos(theta_circle), y = radius_max * sin(theta_circle))

  nose <- data.frame(
    x = c(0, -0.08 * radius_max, 0.08 * radius_max, 0),
    y = c(1.15, 1.02, 1.02, 1.15) * radius_max
  )

  theta_ear <- seq(-pi / 2, pi / 2, length.out = 200)
  ear_left <- data.frame(
    x = -1.05 * radius_max - 0.08 * radius_max * cos(theta_ear),
    y = 0.25 * radius_max * sin(theta_ear)
  )
  ear_right <- data.frame(
    x = 1.05 * radius_max + 0.08 * radius_max * cos(theta_ear),
    y = 0.25 * radius_max * sin(theta_ear)
  )

  clim <- max(abs(topo$diff_change), na.rm = TRUE)
  if (!is.finite(clim) || clim == 0) clim <- 1e-6

  sb <- as.numeric(stats::quantile(topo$abs_change, probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
  sb <- unique(round(sb, 3))
  if (length(sb) < 3) sb <- sort(unique(c(sb, max(topo$abs_change, na.rm = TRUE))))

  ggplot() +
    geom_path(data = circle, aes(x, y), color = "black", linewidth = 1.2) +
    geom_path(data = nose, aes(x, y), color = "black", linewidth = 1.2) +
    geom_path(data = ear_left, aes(x, y), color = "black", linewidth = 1.2) +
    geom_path(data = ear_right, aes(x, y), color = "black", linewidth = 1.2) +
    geom_point(
      data = topo,
      aes(x, y, fill = diff_change, size = abs_change),
      shape = 21, stroke = 0.9, color = "black", alpha = 0.95
    ) +
    geom_point(
      data = subset(topo, sig_fdr),
      aes(x, y),
      shape = 21, size = 12, stroke = 2.0,
      fill = NA, color = "black", alpha = 1
    ) +
    geom_text(
      data = topo,
      aes(x, y, label = label_std),
      vjust = -1.9, size = 3.2, fontface = "bold", color = "black"
    ) +
    scale_fill_gradient2(
      low = "#2B6CB0", mid = "white", high = "#C53030",
      midpoint = 0,
      limits = c(-clim, clim),
      breaks = c(-clim, 0, clim),
      labels = function(x) sprintf("%.2f", x),
      name = "Δ rate\n(N−P)"
    ) +
    scale_size_continuous(
      range = c(3.2, 11),
      breaks = sb,
      labels = function(x) sprintf("%.2f", x),
      name = "|Δ|"
    ) +
    guides(
      fill = guide_colorbar(
        barheight = unit(3.2, "cm"),
        barwidth = unit(0.35, "cm"),
        ticks.colour = "black",
        frame.colour = "black",
        title.position = "top",
        label.position = "right",
        order = 1
      ),
      size = guide_legend(
        override.aes = list(shape = 21, fill = "white", color = "black", alpha = 1),
        title.position = "top",
        order = 2
      )
    ) +
    coord_equal(
      xlim = c(-radius_max * 1.45, radius_max * 1.75),
      ylim = c(-radius_max * 1.25, radius_max * 1.25),
      clip = "off"
    ) +
    labs(
      title = "B. LowAlpha Top-Z hub change topography",
      subtitle = paste0(
        "Outer ring indicates channel-level group difference in change (Welch test), BH-FDR q<",
        alpha_fdr, "."
      )
    ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.15, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 11, hjust = 0.15, margin = margin(b = 10)),
      legend.position = c(0.86, 0.42),
      legend.justification = c(0, 0.5),
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.spacing.y = unit(0.15, "cm"),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.margin = margin(t = 10, r = 55, b = 10, l = 10)
    )
}

panel_a <- create_panel_a(stats_file)
panel_b <- create_panel_b(nodelevel_file, locs_file, top_prop = TOP_PROP)

fig2 <- panel_a + panel_b + patchwork::plot_layout(widths = c(1, 1.25))

png_file <- file.path(output_dir, "Figure2_TopZ_LowAlpha_SpatialLocalization.png")
pdf_file <- file.path(output_dir, "Figure2_TopZ_LowAlpha_SpatialLocalization.pdf")

ggsave(png_file, fig2, width = 14, height = 6, dpi = 600, bg = "white")
ggsave(pdf_file, fig2, width = 14, height = 6, device = cairo_pdf)

cat("Saved:\n", png_file, "\n", pdf_file, "\n")
