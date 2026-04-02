# =============================
# Figure 3B — Path diagram (single mediator, THREE outcomes)
# ONLY CHANGED:
# 1) save to D:/r/nic/3b
# 2) enlarge text inside boxes ONLY
# 3) fix duplicated hjust error
# =============================

library(tidyverse)
library(readxl)
library(janitor)
library(ggtext)
library(grid)

# ---------- PATHS ----------
path_std  <- "D:/r/nic/mediation_results_main_split/main_split__TOPZ_BOTH_Genuine_key_results.csv"
path_dict <- "D:/r/nic修改/csv/NeuroEPO PD Motor & Cognition Variable Dictionary_notitle.xlsx"
out_png   <- "D:/r/nic/3b/Figure3B_SEM_Dose_HighAlphaTemporal_TopZ_ThreePathways_Fixed.png"

# ---------- CREATE OUTPUT DIR ----------
out_dir <- dirname(out_png)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ---------- HELPERS ----------
canon_abbrev <- function(x) {
  x %>%
    sub("_(post|pre)$", "", ., ignore.case = TRUE) %>%
    sub("_a\\d+_\\d+$", "", ., ignore.case = TRUE) %>%
    sub("_a\\d+$", "", ., ignore.case = TRUE)
}

adj_dir <- function(x, lower_better) if (isTRUE(lower_better)) -x else x

short_names <- c(
  "pronesupinelh" = "Pronation–supination\n(Left)",
  "agilityll"     = "Leg agility\n(Left)",
  "taplf"         = "Toe tapping\n(Left)"
)

# ---------- READ DATA ----------
sem_data <- readr::read_csv(path_std, show_col_types = FALSE)

dict <- readxl::read_excel(path_dict, sheet = 1) %>%
  janitor::clean_names() %>%
  transmute(
    abbrev    = str_trim(abbreviation),
    direction = str_trim(direction),
    full_name = str_trim(full_english_name)
  )

get_lower_better <- function(y_post) {
  key <- canon_abbrev(y_post)
  dir <- dict$direction[match(key, dict$abbrev)]
  if (is.na(dir)) stop("Missing direction for: ", y_post)
  str_to_lower(dir) %>% str_detect("^lower is better")
}

get_label <- function(y_post) {
  key <- canon_abbrev(y_post)
  lab <- unname(short_names[key])
  if (!is.na(lab) && lab != "") return(lab)
  full <- dict$full_name[match(key, dict$abbrev)]
  if (!is.na(full) && full != "") return(full)
  key
}

# ---------- SELECT OUTCOMES ----------
y1 <- "pronesupinelh_post"
y2 <- "agilityll_post"
y3 <- "taplf_post"

path1 <- sem_data %>%
  filter(
    mediator_metric == "topz_rate",
    mediator_band == "HighAlpha",
    mediator_region == "Temporal",
    y_post == y1
  )

path2 <- sem_data %>%
  filter(
    mediator_metric == "topz_rate",
    mediator_band == "HighAlpha",
    mediator_region == "Temporal",
    y_post == y2
  )

path3 <- sem_data %>%
  filter(
    mediator_metric == "topz_rate",
    mediator_band == "HighAlpha",
    mediator_region == "Temporal",
    y_post == y3
  )

stopifnot(nrow(path1) == 1, nrow(path2) == 1, nrow(path3) == 1)

lower_better1 <- get_lower_better(path1$y_post)
lower_better2 <- get_lower_better(path2$y_post)
lower_better3 <- get_lower_better(path3$y_post)

label1 <- get_label(path1$y_post)
label2 <- get_label(path2$y_post)
label3 <- get_label(path3$y_post)

a <- mean(c(path1$a_est_std, path2$a_est_std, path3$a_est_std), na.rm = TRUE)

# ---------- NODES ----------
mediator_label <- "HighAlpha Temporal\nTop\u2011Z rate"

nodes <- tibble(
  node = c("Dose", mediator_label, label1, label2, label3),
  x    = c(0.00, 1.00, 2.10, 2.10, 2.10),
  y    = c(0.80, 0.80, 1.20, 0.80, 0.40),
  type = c("Dose", "Mediator", "Outcome", "Outcome", "Outcome")
)

node_x <- function(nm) nodes$x[match(nm, nodes$node)]
node_y <- function(nm) nodes$y[match(nm, nodes$node)]

# ---------- EDGES ----------
edges <- tibble(
  from = c("Dose", mediator_label, mediator_label, mediator_label),
  to   = c(mediator_label, label1, label2, label3),
  
  # 原始节点坐标
  x_from = node_x(from),
  y_from = node_y(from),
  x_to   = node_x(to),
  y_to   = node_y(to),
  
  # 调整箭头起点：从 mediator 方框的边缘开始
  xstart = if_else(from == "Dose", x_from, x_from + 0.35),
  ystart = case_when(
    to == label1 ~ y_from + 0.10,
    to == label2 ~ y_from,
    to == label3 ~ y_from - 0.10,
    TRUE ~ y_from
  ),
  
  # 调整箭头终点：到 outcome 左边缘
  xend = if_else(from == "Dose", x_to - 0.35, x_to - 0.35),
  yend = y_to,
  
  curvature = c(0.00, 0.15, 0.00, -0.15),
  
  label = c(
    paste0("a = ", sprintf("%.2f", a)),
    paste0("b\u2081 = ", sprintf("%.2f", path1$b_est_std),
           "\nind\u2081 = ", sprintf("%.2f", adj_dir(path1$ind_est_std, lower_better1))),
    paste0("b\u2082 = ", sprintf("%.2f", path2$b_est_std),
           "\nind\u2082 = ", sprintf("%.2f", adj_dir(path2$ind_est_std, lower_better2))),
    paste0("b\u2083 = ", sprintf("%.2f", path3$b_est_std),
           "\nind\u2083 = ", sprintf("%.2f", adj_dir(path3$ind_est_std, lower_better3)))
  ),
  
  # 标签位置（保持原始）
  lx = c(0.48, 1.58, 1.58, 1.58),
  ly = c(0.84, 1.08, 0.80, 0.56)
)

# ---------- PLOT ----------
p_path <- ggplot() +
  # 保持你的原始连线写法
  geom_curve(
    data = edges,
    aes(
      x = xstart, y = ystart,
      xend = xend, yend = yend,
      curvature = curvature
    ),
    arrow = arrow(length = unit(0.28, "cm"), type = "closed"),
    linewidth = 1.2,
    color = "gray25"
  ) +
  
  # 边标签保持原样
  geom_text(
    data = edges,
    aes(x = lx, y = ly, label = label),
    size = 4.2,
    fontface = "bold",
    lineheight = 0.9
  ) +
  
  # 节点方框：完全保持原始大小
  geom_rect(
    data = nodes,
    aes(
      xmin = x - 0.35, xmax = x + 0.35,
      ymin = y - 0.12, ymax = y + 0.12,
      fill = type
    ),
    color = "black",
    linewidth = 0.8
  ) +
  
  # 只把方框里的字调大
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = node),
    size = 6.8,
    fontface = "bold",
    lineheight = 0.85
  ) +
  
  scale_fill_manual(
    values = c("Dose" = "#C9CED3", "Mediator" = "#F1C40F", "Outcome" = "#2E86C1"),
    guide = "none"
  ) +
  
  coord_cartesian(xlim = c(-0.30, 2.55), ylim = c(0.20, 1.40)) +
  
  labs(
    title = "<b>Complete High\u2011Alpha Temporal Mediation Pathways</b>",
    subtitle = "All three dose–outcome pathways mediated by high-alpha temporal top-Z rate reorganization (standardized coefficients)"
  ) +
  
  theme_void() +
  theme(
    plot.margin = margin(t = 10, r = 8, b = 6, l = 8),
    plot.title = ggtext::element_markdown(hjust = 0.5, size = 16.5),
    plot.subtitle = ggtext::element_markdown(hjust = 0.5, size = 11)
  )

ggsave(out_png, p_path, width = 14, height = 6.5, dpi = 600, bg = "white")

cat("✓ Figure saved:", out_png, "\n")

p_path