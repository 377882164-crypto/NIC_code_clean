suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(readr)
  library(igraph)
  library(broom.mixed)
  library(lme4)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(rlang)
  glmm_ok <- requireNamespace("glmmTMB", quietly = TRUE)
})

# ============================================================
# NIC clean mainline only: EEG network outputs -> TopZ -> GLMM -> TOPZ figure
# Built by trimming the original script WITHOUT changing:
#   1) TopZ definition
#   2) TopZ rate logic
#   3) GLMM formula / factor reference / data filtering
# ============================================================

# ====== paths ======
out_root <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw/05_results/clean_topz_integrated"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(out_root, "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
topz_dir <- file.path(out_root, "topz_results")
dir.create(topz_dir, showWarnings = FALSE, recursive = TRUE)
topz_fig_dir <- file.path(fig_dir, "topz")
dir.create(topz_fig_dir, showWarnings = FALSE, recursive = TRUE)

base_dir_placebo_pre  <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw/05_results/network/Placebo_pre"
base_dir_placebo_post <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw/05_results/network/Placebo_post"
base_dir_neuro_pre    <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw/05_results/network/NeuroEPO_pre"
base_dir_neuro_post   <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw/05_results/network/NeuroEPO_post"

file_placebo <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw/01_raw_data/wide_dose_Hz_Placebo.csv"
file_neuro   <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw/01_raw_data/wide_dose_Hz_Neuro.csv"

file_placebo_pre_cen  <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw/02_auxiliary/Multi_Centrality_Results2/placebo_pre_all_centralities.csv"
file_placebo_post_cen <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw/02_auxiliary/Multi_Centrality_Results2/placebo_post_all_centralities.csv"
file_neuro_pre_cen    <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw/02_auxiliary/Multi_Centrality_Results2/neuro_pre_all_centralities.csv"
file_neuro_post_cen   <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw/02_auxiliary/Multi_Centrality_Results2/neuro_post_all_centralities.csv"

# Optional locked old results for strict equality check after clean run.
# Leave as NA_character_ to skip comparison.
old_topz_dir <- NA_character_

# ====== constants copied from old script ======
breaks_fg <- c(0.5, 4, 8, 10, 13, 30)
labels_fg <- c("Delta 0.5--4","Theta 4--8","Low-Alpha 8--10","High-Alpha 10--13","Beta 13--30")
bands_order <- c("Delta","Theta","LowAlpha","HighAlpha","Beta")
regions_order <- c("Frontal","Central","Parietal","Occipital","Temporal")

topz_top_prop <- 0.10
topz_min_nodes <- 5

# ====== helpers copied from old mainline ======
theme_pub <- function(base_size = 12, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = base_size * 0.9),
      axis.text.y = ggplot2::element_text(size = base_size * 0.9),
      strip.text  = ggplot2::element_text(face = "bold"),
      plot.title  = ggplot2::element_text(face = "bold", hjust = 0.0),
      legend.title= ggplot2::element_text(face = "bold")
    )
}

freq_group_lab <- function(freq_num) {
  cut(freq_num, breaks = breaks_fg,
      labels = labels_fg[1:(length(breaks_fg)-1)],
      include.lowest = TRUE, right = FALSE)
}

map_region <- function(ch) {
  s <- as.character(ch)
  case_when(
    str_starts(s, regex("Fp|Fz|F", ignore_case = TRUE)) ~ "Frontal",
    str_starts(s, regex("Cz|C",   ignore_case = TRUE)) ~ "Central",
    str_starts(s, regex("P",      ignore_case = TRUE)) ~ "Parietal",
    str_starts(s, regex("O",      ignore_case = TRUE)) ~ "Occipital",
    str_starts(s, regex("T",      ignore_case = TRUE)) ~ "Temporal",
    TRUE ~ "Other"
  )
}

fmt2 <- function(x) sprintf("%.3f", as.numeric(x))

parse_col <- function(x) {
  s <- sub("^X+", "", x)
  time <- ifelse(str_detect(s, regex("post$", ignore_case = TRUE)), "post",
                 ifelse(str_detect(s, regex("pre$",  ignore_case = TRUE)), "pre", NA_character_))
  s2 <- sub("(?i)(?:_)?(?:pre|post)$","", s, perl = TRUE)
  m <- regexpr("(?i)\\d+(?:\\.\\d+)?(?=\\s*Hz)", s2, perl = TRUE)
  if (m[1] > 0) {
    freq_val <- substr(s2, m[1], m[1] + attr(m,"match.length") - 1)
    channel_raw <- gsub("(?i)\\d+(?:\\.\\d+)?\\s*Hz", "", s2, perl = TRUE)
  } else {
    m2 <- regexpr("\\d+(?:\\.\\d+)?", s2, perl = TRUE)
    if (m2[1] > 0) {
      freq_val <- substr(s2, m2[1], m2[1] + attr(m2,"match.length") - 1)
      channel_raw <- paste0(substr(s2, 1, m2[1]-1), substr(s2, m2[1] + attr(m2, "match.length"), nchar(s2)))
    } else { freq_val <- NA_character_; channel_raw <- s2 }
  }
  channel <- gsub("[^A-Za-z0-9]", "", channel_raw)
  channel <- sub("^[0-9]+", "", channel)
  tibble(orig = x, Time = time, Freq = suppressWarnings(as.numeric(freq_val)), Channel = channel)
}

find_comm_detail_file <- function(dir_path) {
  if (!dir.exists(dir_path)) return(NA_character_)
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(NA_character_)
  pat1 <- "^4\\.1Community_Variables_.*\\.csv$"
  pat2 <- "^Community_Variables_.*\\.csv$"
  cand1 <- files[grepl(pat1, basename(files), ignore.case = FALSE)]
  cand2 <- files[grepl(pat2, basename(files), ignore.case = FALSE)]
  candidates <- if (length(cand1) > 0) cand1 else cand2
  if (length(candidates) == 0) return(NA_character_)
  info <- file.info(candidates)
  candidates[which.max(info$mtime)]
}

normalize_colnames <- function(df) {
  nm <- names(df)
  nm <- gsub("^\ufeff", "", nm)
  nm <- gsub("\\s+", "", nm, perl = TRUE)
  low <- tolower(nm)
  map <- c(
    "communityid"    = "CommunityID",
    "communityindex" = "CommunityIndex",
    "node"           = "Node",
    "channel"        = "Channel",
    "frequency"      = "Frequency",
    "freq"           = "Frequency",
    "freqgroup"      = "FreqGroup",
    "strength"       = "Strength",
    "betweenness"    = "Betweenness",
    "region"         = "Region"
  )
  std <- ifelse(low %in% names(map), map[low], nm)
  names(df) <- std
  df
}

try_fetch_centrality <- function(Group_label, Time_label) {
  cen_file <- switch(
    paste0(Group_label, "_", Time_label),
    "Placebo_pre"  = file_placebo_pre_cen,
    "Placebo_post" = file_placebo_post_cen,
    "NeuroEPO_pre" = file_neuro_pre_cen,
    "NeuroEPO_post"= file_neuro_post_cen,
    NULL
  )
  if (is.null(cen_file) || !file.exists(cen_file)) {
    warning("Centrality file not found for ", Group_label, " ", Time_label, ": ", cen_file)
    return(NULL)
  }
  cen <- suppressMessages(readr::read_csv(cen_file, show_col_types = FALSE))
  cen <- normalize_colnames(cen)
  if (!("Node" %in% names(cen))) {
    if ("node" %in% names(cen)) {
      names(cen)[names(cen) == "node"] <- "Node"
    } else {
      warning("centrality file 缺少 Node 列: ", cen_file)
      return(NULL)
    }
  }
  strength_col <- names(cen)[names(cen) %in% c("Strength","Strength_mm","Strength_z")]
  betw_col     <- names(cen)[names(cen) %in% c("Betweenness","Betweenness_z")]
  degree_col   <- names(cen)[names(cen) %in% c("Degree","degree","Degree_z")]
  strength_col <- strength_col[1]
  betw_col     <- betw_col[1]
  degree_col   <- degree_col[1]
  if (is.na(strength_col) && is.na(betw_col) && is.na(degree_col)) {
    warning("centrality file 中没有 Strength/Betweenness/Degree 相关列: ", cen_file)
    return(NULL)
  }
  cen %>%
    transmute(
      Node = .data[["Node"]],
      Strength    = if (!is.na(strength_col)) as.numeric(.data[[strength_col]]) else NA_real_,
      Betweenness = if (!is.na(betw_col))     as.numeric(.data[[betw_col]])     else NA_real_,
      Degree      = if (!is.na(degree_col))   as.numeric(.data[[degree_col]])   else NA_real_
    )
}

identify_topz_nodes_per_person <- function(df_person, id, group, time) {
  dfp_all <- df_person %>% filter(is.finite(Value_z))
  out_list <- list()
  get_top_idx <- function(n) max(1, ceiling(n * topz_top_prop))
  if (nrow(dfp_all) >= topz_min_nodes) {
    k_all <- get_top_idx(nrow(dfp_all))
    idx_all <- order(abs(dfp_all$Value_z), decreasing = TRUE)
    top_idx_all <- idx_all[seq_len(k_all)]
    out_list[["All"]] <- dfp_all[top_idx_all, ] %>% mutate(TopZScope = "All")
  }
  by_band <- df_person %>% filter(!is.na(Band), is.finite(Value_z)) %>% group_split(Band, .keep = TRUE)
  names(by_band) <- df_person %>% filter(!is.na(Band), is.finite(Value_z)) %>% pull(Band) %>% unique()
  for (b in names(by_band)) {
    db <- by_band[[b]]
    if (nrow(db) >= topz_min_nodes) {
      k_b <- get_top_idx(nrow(db))
      idx_b <- order(abs(db$Value_z), decreasing = TRUE)
      top_idx_b <- idx_b[seq_len(k_b)]
      out_list[[paste0("Band_", b)]] <- db[top_idx_b, ] %>% mutate(TopZScope = paste0("Band_", b))
    }
  }
  if (length(out_list) == 0) return(tibble())
  bind_rows(out_list) %>%
    mutate(is_topz = 1L) %>%
    select(TopZScope, Node, Channel, Frequency, Band, Region, Value_z, is_topz)
}

summarise_topz_rates_per_person <- function(df_person, topz_person) {
  base_tab <- df_person %>%
    filter(!is.na(Band), !is.na(Region), is.finite(Value_z)) %>%
    group_by(Band, Region) %>%
    summarise(n_nodes = n_distinct(Node), .groups = "drop")
  if (nrow(base_tab) == 0) return(tibble())
  topz_tab <- topz_person %>%
    filter(startsWith(TopZScope, "Band_")) %>%
    mutate(Band = sub("^Band_", "", TopZScope)) %>%
    group_by(Band, Region) %>%
    summarise(n_topz = n_distinct(Node), .groups = "drop")
  base_tab %>%
    left_join(topz_tab, by = c("Band", "Region")) %>%
    mutate(n_topz = replace_na(n_topz, 0L),
           topz_rate = if_else(n_nodes > 0, n_topz / n_nodes, NA_real_))
}

fit_glmm_binom <- function(formula, data) {
  if (!glmm_ok) return(list(model = NULL, family = "none"))
  fit1 <- tryCatch(glmmTMB::glmmTMB(formula, data = data, family = binomial(), REML = FALSE),
                   error = function(e) NULL, warning = function(w) NULL)
  if (!is.null(fit1)) return(list(model = fit1, family = "binomial"))
  fit2 <- tryCatch(glmmTMB::glmmTMB(formula, data = data,
                                    family = glmmTMB::beta_family(link = "logit"),
                                    REML = FALSE),
                   error = function(e) NULL, warning = function(w) NULL)
  if (!is.null(fit2)) return(list(model = fit2, family = "beta"))
  list(model = NULL, family = "none")
}

tidy_glmm <- function(fit_obj) {
  if (is.null(fit_obj$model)) return(NULL)
  tab <- tryCatch(broom.mixed::tidy(fit_obj$model, effects = "fixed", component = "cond"),
                  error = function(e) NULL)
  if (is.null(tab)) return(NULL)
  if (!"p.value" %in% names(tab)) tab$p.value <- NA_real_
  tab$family_used <- fit_obj$family
  tab
}

make_individual_metrics <- function(base_dir, file_raw, Group_label, Time_label, out_root) {
  out_dir <- file.path(out_root, paste0(Group_label, "_", Time_label))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  file_comm_detail <- find_comm_detail_file(base_dir)
  if (is.na(file_comm_detail) || !file.exists(file_comm_detail)) {
    warning("Community detail CSV not found: ", base_dir)
    return(list(topz_nodes = tibble(), topz_rates = tibble(), detail_with_cent = tibble(), out_dir = out_dir))
  }
  message("Using community detail file: ", file_comm_detail)
  comm_detail <- suppressMessages(readr::read_csv(file_comm_detail, show_col_types = FALSE))
  comm_detail <- normalize_colnames(comm_detail)
  if (!all(c("Node","CommunityIndex") %in% names(comm_detail))) {
    stop("Community detail lacks required columns (Node or CommunityIndex). Existing: ",
         paste(names(comm_detail), collapse = ", "))
  }

  sep <- tidyr::separate_wider_delim(tibble(Node = comm_detail$Node),
                                     Node, delim = "_",
                                     names = c("Channel_tmp", "Frequency_str"), too_few = "align_start")
  comm_detail$Channel   <- sep$Channel_tmp
  comm_detail$Frequency <- suppressWarnings(as.numeric(sep$Frequency_str))
  comm_detail$FreqGroup <- as.character(freq_group_lab(comm_detail$Frequency))
  comm_detail$Region    <- map_region(comm_detail$Channel)

  if (!("Strength" %in% names(comm_detail)))    comm_detail$Strength <- NA_real_
  if (!("Betweenness" %in% names(comm_detail))) comm_detail$Betweenness <- NA_real_
  if (!("Degree" %in% names(comm_detail)))      comm_detail$Degree <- NA_real_

  cen_fill <- try_fetch_centrality(Group_label, Time_label)
  if (!is.null(cen_fill)) {
    comm_detail <- comm_detail %>%
      left_join(cen_fill, by = "Node", suffix = c("", "_fromCen")) %>%
      mutate(
        Strength    = ifelse(is.na(Strength)    & !is.na(Strength_fromCen),    Strength_fromCen,    Strength),
        Betweenness = ifelse(is.na(Betweenness) & !is.na(Betweenness_fromCen), Betweenness_fromCen, Betweenness),
        Degree      = ifelse(is.na(Degree)      & !is.na(Degree_fromCen),      Degree_fromCen,      Degree)
      ) %>%
      select(-ends_with("_fromCen"))
  }

  raw <- suppressMessages(readr::read_csv(file_raw, show_col_types = FALSE))
  nm0 <- names(raw)
  nm_clean <- gsub("^\ufeff", "", nm0)
  nm_clean <- gsub("\\s+", "", nm_clean, perl = TRUE)
  nm_lower <- tolower(nm_clean)
  std_names <- nm0
  std_names[nm_lower %in% c("group","grp","treatment","arm")] <- "GroupCode"
  std_names[nm_lower %in% c("id","subject","subj","participant")] <- "ID"
  std_names[nm_lower %in% c("dose","dose_level","dosage")] <- "dose"
  names(raw) <- std_names
  if (any(duplicated(names(raw)))) names(raw) <- make.unique(names(raw))

  key_cols <- intersect(c("ID","dose","GroupCode"), names(raw))
  if (!("ID" %in% key_cols)) stop("Raw wide table lacks required column ID. Current names: ", paste(names(raw), collapse = ", "))
  meta <- purrr::map_dfr(names(raw), parse_col)

  node_info <- comm_detail %>%
    transmute(Node, Channel, Frequency, FreqGroup, Region, CommunityIndex,
              NodeKey = paste0(Channel, "_", fmt2(Frequency))) %>%
    normalize_colnames()

  meta_time <- meta %>%
    filter(Time == !!Time_label & !is.na(Freq) & !is.na(Channel)) %>%
    mutate(NodeKey = paste0(Channel, "_", fmt2(Freq)))

  col_to_node_min <- meta_time %>%
    inner_join(node_info %>% select(NodeKey, Node, CommunityIndex), by = "NodeKey") %>%
    select(orig, Node, CommunityIndex)

  if (nrow(col_to_node_min) == 0) {
    stop(Group_label, " ", Time_label, ": cannot map raw columns to nodes (NodeKey mismatch).")
  }

  col_to_node <- col_to_node_min %>%
    left_join(node_info %>% select(Node, Channel, Frequency, FreqGroup, Region, CommunityIndex),
              by = c("Node","CommunityIndex"))

  dat_time <- raw %>% select(all_of(key_cols), all_of(unique(col_to_node$orig)))

  long_ind <- dat_time %>%
    pivot_longer(cols = -all_of(key_cols), names_to = "orig", values_to = "Value") %>%
    left_join(col_to_node, by = "orig") %>%
    filter(!is.na(Node), !is.na(CommunityIndex)) %>%
    group_by(ID) %>%
    mutate(Value_z = as.numeric(scale(Value))) %>%
    ungroup() %>%
    mutate(Group = Group_label, Time = Time_label)

  detail_with_cent <- long_ind %>%
    mutate(
      Strength_group    = comm_detail$Strength[match(Node, comm_detail$Node)],
      Betweenness_group = comm_detail$Betweenness[match(Node, comm_detail$Node)],
      Degree_group      = comm_detail$Degree[match(Node, comm_detail$Node)],
      Band = case_when(
        str_detect(FreqGroup, regex("Delta", ignore_case = TRUE)) ~ "Delta",
        str_detect(FreqGroup, regex("Theta", ignore_case = TRUE)) ~ "Theta",
        str_detect(FreqGroup, regex("Low-Alpha", ignore_case = TRUE)) ~ "LowAlpha",
        str_detect(FreqGroup, regex("High-Alpha", ignore_case = TRUE)) ~ "HighAlpha",
        str_detect(FreqGroup, regex("Beta", ignore_case = TRUE)) ~ "Beta",
        TRUE ~ NA_character_
      )
    )

  person_nodes <- detail_with_cent %>%
    select(Group, ID, Time, Node, Channel, Frequency, Band, Region, Value_z) %>%
    distinct()

  person_keys <- person_nodes %>%
    distinct(Group, ID, Time)

  person_groups <- person_nodes %>%
    group_by(Group, ID, Time) %>%
    group_split(.keep = TRUE)

  topz_nodes_person <- purrr::map2_dfr(person_groups, seq_along(person_groups), function(dfp, i) {
    key <- person_keys[i, , drop = FALSE]
    grp <- as.character(key$Group[[1]])
    person_id <- key$ID[[1]]
    tm <- as.character(key$Time[[1]])

    res <- identify_topz_nodes_per_person(dfp, id = person_id, group = grp, time = tm)
    if (nrow(res) == 0) {
      return(tibble())
    }

    res %>%
      mutate(Group = grp, ID = person_id, Time = tm, .before = 1)
  })

  topz_rates_person <- purrr::map2_dfr(person_groups, seq_along(person_groups), function(dfp, i) {
    key <- person_keys[i, , drop = FALSE]
    grp <- as.character(key$Group[[1]])
    person_id <- key$ID[[1]]
    tm <- as.character(key$Time[[1]])

    tz_p <- topz_nodes_person %>%
      filter(Group == grp, ID == person_id, Time == tm)

    res <- summarise_topz_rates_per_person(dfp, tz_p)
    if (nrow(res) == 0) {
      return(tibble())
    }

    res %>%
      mutate(Group = grp, ID = person_id, Time = tm, .before = 1)
  })

  write_csv(topz_nodes_person, file.path(out_dir, "topz_nodes_person.csv"))
  write_csv(topz_rates_person, file.path(out_dir, "topz_rates_person.csv"))
  write_csv(detail_with_cent %>%
              select(Group, Time, ID, Node, Channel, Frequency, Band, Region,
                     Value_z, Strength_group, Betweenness_group, Degree_group),
            file.path(out_dir, "signal_activity_nodelevel_all.csv"))

  list(topz_nodes = topz_nodes_person, topz_rates = topz_rates_person, detail_with_cent = detail_with_cent, out_dir = out_dir)
}

# ====== run four batches ======
batches <- tribble(
  ~Group_label, ~Time_label, ~base_dir,             ~file_raw,
  "Placebo",    "pre",       base_dir_placebo_pre,  file_placebo,
  "Placebo",    "post",      base_dir_placebo_post, file_placebo,
  "NeuroEPO",   "pre",       base_dir_neuro_pre,    file_neuro,
  "NeuroEPO",   "post",      base_dir_neuro_post,   file_neuro
)

results <- pmap(batches, function(Group_label, Time_label, base_dir, file_raw) {
  make_individual_metrics(base_dir, file_raw, Group_label, Time_label, out_root)
})

# ====== aggregate exact TopZ outputs ======
topz_nodes_all <- tryCatch(map_dfr(results, "topz_nodes"), error = function(e) tibble())
topz_rates_all <- tryCatch(map_dfr(results, "topz_rates"), error = function(e) tibble())
detail_all     <- tryCatch(map_dfr(results, "detail_with_cent"), error = function(e) tibble())

if (nrow(topz_nodes_all) > 0) {
  topz_nodes_all <- topz_nodes_all %>%
    mutate(
      Band   = if ("Band" %in% names(.)) factor(as.character(Band), levels = bands_order) else NA,
      Region = if ("Region" %in% names(.)) factor(as.character(Region), levels = regions_order) else NA,
      Group  = factor(as.character(Group), levels = c("Placebo","NeuroEPO")),
      Time   = factor(as.character(Time),  levels = c("pre","post"))
    )
  write_csv(topz_nodes_all, file.path(topz_dir, "topz_nodes_all.csv"))
}

if (nrow(topz_rates_all) > 0) {
  topz_rates_all <- topz_rates_all %>%
    mutate(
      Band   = factor(as.character(Band),   levels = bands_order),
      Region = factor(as.character(Region), levels = regions_order),
      Group  = factor(as.character(Group),  levels = c("Placebo","NeuroEPO")),
      Time   = factor(as.character(Time),   levels = c("pre","post"))
    )
  write_csv(topz_rates_all, file.path(topz_dir, "topz_rates_all.csv"))
}

if (nrow(detail_all) > 0) {
  detail_all <- detail_all %>%
    mutate(
      Group = factor(as.character(Group), levels = c("Placebo","NeuroEPO")),
      Time  = factor(as.character(Time),  levels = c("pre","post")),
      Band  = factor(as.character(Band),  levels = bands_order),
      Region= factor(as.character(Region),levels = regions_order)
    )
  write_csv(
    detail_all %>% select(Group, Time, ID, Node, Channel, Frequency, Band, Region,
                          Value_z, Strength_group, Betweenness_group, Degree_group),
    file.path(out_root, "signal_activity_nodelevel_all.csv")
  )
}

# ====== GLMM copied from old TopZ branch ======
if (nrow(topz_rates_all) > 0) {
  do_topz_glmm <- function(df, band, region) {
    sub <- df %>% filter(as.character(Band) == band, as.character(Region) == region)
    if (nrow(sub) == 0 || n_distinct(sub$ID) < 3) return(NULL)
    dd <- sub %>%
      mutate(
        k = as.integer(round(n_topz)),
        n = as.integer(round(n_nodes))
      ) %>%
      filter(is.finite(k), is.finite(n), n > 0, k >= 0, k <= n)
    if (nrow(dd) < 8) return(NULL)
    f <- cbind(k, n - k) ~ Group*Time + (1|ID)
    fito <- fit_glmm_binom(f, dd)
    if (is.null(fito$model)) {
      y0 <- qlogis(pmin(pmax(dd$k/dd$n, 1e-6), 1-1e-6))
      fit_lmm <- tryCatch(lme4::lmer(y0 ~ Group*Time + (1|ID), data = dd,
                                     control = lme4::lmerControl(optimizer = "bobyqa")),
                          error = function(e) NULL)
      if (is.null(fit_lmm)) return(NULL)
      tab <- broom.mixed::tidy(fit_lmm, effects = "fixed")
      if (!"p.value" %in% names(tab)) tab$p.value <- NA_real_
      tab$family_used <- "LMM-logit-fallback"
      return(tab %>% mutate(Band = band, Region = region))
    }
    tab <- tidy_glmm(fito)
    if (is.null(tab)) return(NULL)
    tab %>% mutate(Band = band, Region = region)
  }

  combos <- tidyr::expand_grid(Band = bands_order, Region = regions_order)
  topz_glmm_tabs <- purrr::pmap_dfr(list(as.character(combos$Band), as.character(combos$Region)),
                                    ~do_topz_glmm(topz_rates_all, ..1, ..2))
  if (nrow(topz_glmm_tabs) > 0) {
    write_csv(topz_glmm_tabs, file.path(topz_dir, "stats_GLMM_topz_all_terms.csv"))
    topz_glmm_int <- topz_glmm_tabs %>%
      filter(term == "GroupNeuroEPO:Timepost") %>%
      group_by(Band, Region) %>%
      mutate(p_adj = ifelse(is.finite(p.value), p.adjust(p.value, "BH"), NA_real_)) %>%
      ungroup()
    write_csv(topz_glmm_int, file.path(topz_dir, "stats_GLMM_topz_interaction.csv"))

    if (nrow(topz_glmm_int) > 0) {
      show_tab <- topz_glmm_int %>%
        mutate(Band = factor(as.character(Band), levels = bands_order),
               Region = factor(as.character(Region), levels = regions_order))
      p_tz <- ggplot(show_tab, aes(Region, Band, fill = -log10(p_adj))) +
        geom_tile(color = "grey90") +
        scale_fill_gradient(low = "#eef5ff", high = "#084594", na.value = "grey95") +
        geom_text(aes(label = ifelse(is.finite(p_adj) & p_adj < 0.05, "*", "")), size = 5) +
        labs(title = "Interaction (Group × Time) on Top‑Z rate (GLMM)", fill = "-log10(FDR)") +
        theme_pub(12) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
      ggsave(file.path(topz_fig_dir, "topzrate_crossgroup_interaction_heatmap_GLMM.png"),
             p_tz, width = 10.5, height = 7, dpi = 360, bg = "white")
      ggsave(file.path(out_root, "topzrate_crossgroup_interaction_heatmap_GLMM.png"),
             p_tz, width = 10.5, height = 7, dpi = 360, bg = "white")
    }
  } else {
    warning("top‑Z GLMM: 样本或方差不足，或模型未能收敛。")
  }
}

# ====== optional strict equality check against locked old results ======
compare_exact_csv <- function(clean_file, old_file, key_cols, numeric_cols = NULL) {
  if (!file.exists(clean_file) || !file.exists(old_file)) return(tibble(ok = FALSE, reason = "missing_file"))
  clean <- suppressMessages(readr::read_csv(clean_file, show_col_types = FALSE))
  old   <- suppressMessages(readr::read_csv(old_file, show_col_types = FALSE))
  clean <- clean %>% mutate(across(any_of(key_cols), as.character)) %>% arrange(across(any_of(key_cols)))
  old   <- old   %>% mutate(across(any_of(key_cols), as.character)) %>% arrange(across(any_of(key_cols)))
  if (!identical(names(clean), names(old))) return(tibble(ok = FALSE, reason = "different_columns"))
  if (!identical(clean[key_cols], old[key_cols])) return(tibble(ok = FALSE, reason = "different_keys"))
  if (is.null(numeric_cols)) numeric_cols <- intersect(names(clean), names(old))
  numeric_cols <- numeric_cols[sapply(clean[numeric_cols], is.numeric) & sapply(old[numeric_cols], is.numeric)]
  non_num_cols <- setdiff(names(clean), numeric_cols)
  same_non_num <- identical(clean[non_num_cols], old[non_num_cols])
  max_diff <- 0
  if (length(numeric_cols) > 0) {
    max_diff <- max(abs(as.matrix(clean[numeric_cols]) - as.matrix(old[numeric_cols])), na.rm = TRUE)
    if (!is.finite(max_diff)) max_diff <- 0
  }
  tibble(ok = isTRUE(same_non_num) && max_diff == 0, reason = ifelse(isTRUE(same_non_num) && max_diff == 0, "exact_match", "value_difference"), max_diff = max_diff)
}

if (is.character(old_topz_dir) && length(old_topz_dir) == 1 && !is.na(old_topz_dir)) {
  cmp_topz <- compare_exact_csv(
    clean_file = file.path(topz_dir, "topz_rates_all.csv"),
    old_file   = file.path(old_topz_dir, "topz_rates_all.csv"),
    key_cols   = c("Group", "ID", "Time", "Band", "Region")
  )
  write_csv(cmp_topz, file.path(topz_dir, "clean_vs_old_topz_rates_check.csv"))
  if (!isTRUE(cmp_topz$ok[1])) {
    warning("clean vs old topz_rates_all.csv NOT exact match. Check: ", file.path(topz_dir, "clean_vs_old_topz_rates_check.csv"))
  } else {
    message("✅ clean vs old topz_rates_all.csv exact match")
  }
}

message("Completed: clean TOPZ mainline only. Outputs in: ", normalizePath(out_root))
