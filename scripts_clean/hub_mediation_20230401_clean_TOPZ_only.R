suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(furrr)
  library(progressr)
})

# ===============================================================
# NIC mediation (clean, TOPZ-only)
# Based on the original hub_mediation_20230401.R
# Mainline preserved:
#   topz_rates_all.csv + wide_dose_all3_cleaned_noHz.csv
#   -> Band x Region topz mediation tables
#   -> split-sample bootstrap outputs
# Classic hub branch removed intentionally.
# ===============================================================

# ---------------- path helpers ----------------
project_root <- "C:/Users/11227/OneDrive/Desktop/NIC_reproduction_from_raw"

pick_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit)) {
    stop(
      "None of the candidate files exist. Checked:\n",
      paste(paths, collapse = "\n"),
      call. = FALSE
    )
  }
  hit
}

f_topz <- pick_first_existing(c(
  file.path(project_root, "05_results", "clean_topz_integrated", "topz_results", "topz_rates_all.csv"),
  file.path(project_root, "05_results", "old_topz_integrated",   "topz_results", "topz_rates_all.csv"),
  file.path(project_root, "05_results", "topz_results",          "topz_rates_all.csv"),
  file.path(project_root, "data_derived",                         "topz_rates_all.csv")
))

f_wide <- pick_first_existing(c(
  file.path(project_root, "01_raw_data", "wide_dose_all3_cleaned_noHz.csv")
))

# ---------------- output ----------------
out_dir <- file.path(project_root, "05_results", "mediation_results_main_split_clean_topz")
file_prefix <- "main_split__"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Using topz file: ", f_topz)
message("Using wide file: ", f_wide)
message("Output dir: ", out_dir)

# ---------------- read data ----------------
topz_all <- readr::read_csv(f_topz, show_col_types = FALSE) %>% mutate(ID = as.character(ID))
wide     <- readr::read_csv(f_wide, show_col_types = FALSE) %>% mutate(ID = as.character(ID))

# ---------------- mediation config ----------------
use_mediators <- c("topz_rate")

bands   <- c("Delta", "Theta", "LowAlpha", "HighAlpha", "Beta")
regions <- c("Frontal", "Central", "Parietal", "Occipital", "Temporal")

# ---------------- build mediator wide table ----------------
make_med_table <- function(df, value_col) {
  df %>%
    filter(Band %in% bands, Region %in% regions, Time %in% c("pre", "post")) %>%
    select(ID, Time, Band, Region, !!sym(value_col)) %>%
    mutate(var = paste0(value_col, "__", Band, "__", Region)) %>%
    select(ID, Time, var, value = !!sym(value_col)) %>%
    pivot_wider(names_from = c(Time, var), values_from = value) %>%
    relocate(ID)
}

med_topz <- make_med_table(topz_all, "topz_rate")

# merge to wide
# NOTE: original script used full_join across hub/topz branches.
# Since only TOPZ is retained, we merge only TOPZ branch into wide.
dat0 <- wide %>% left_join(med_topz, by = "ID")

# ---------------- identify Y(post) and paired pre ----------------
y_post_all <- names(dat0) %>%
  keep(~ str_ends(.x, "_post")) %>%
  setdiff(c("age_post", "education_post"))

y_map <- tibble(
  y_post = y_post_all,
  y_pre  = sub("_post$", "_pre", y_post_all)
) %>%
  filter(y_pre %in% names(dat0))

# ---------------- mediator definitions ----------------
mk_mediator_defs <- function(value_col) {
  expand_grid(Band = bands, Region = regions) %>%
    transmute(
      mediator_metric = value_col,
      m_pre  = paste0("pre_",  value_col, "__", Band, "__", Region),
      m_post = paste0("post_", value_col, "__", Band, "__", Region),
      mediator_band = Band,
      mediator_region = Region
    )
}

mediators_def <- map_df(use_mediators, mk_mediator_defs) %>%
  filter(m_pre %in% names(dat0), m_post %in% names(dat0))

# ---------------- utilities ----------------
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

tidy_augment <- function(fit, conf_level = 0.95) {
  ti <- broom::tidy(fit, conf.int = TRUE, conf.level = conf_level)
  gl <- broom::glance(fit)
  list(tidy = ti, glance = gl)
}

zscore_df <- function(df, cols) {
  out <- df
  for (nm in cols) {
    v <- df[[nm]]
    if (is.numeric(v)) {
      s <- stats::sd(v, na.rm = TRUE)
      if (is.finite(s) && s > 0) {
        out[[nm]] <- (v - mean(v, na.rm = TRUE)) / s
      } else {
        out[[nm]] <- 0
      }
    }
  }
  out
}

# Bootstrap: M and Y sampled separately; M model does not require outcome completeness.
boot_indirect_split <- function(dfM, dfY, y_post, y_pre, m_post, m_pre, R = 2000, seed = 2025) {
  set.seed(seed)

  nM <- nrow(dfM)
  nY <- nrow(dfY)

  # unstandardized fits
  fit_m <- lm(reformulate(c("dose", m_pre, "age_post", "education_post"), response = m_post), data = dfM)
  fit_y <- lm(reformulate(c("dose", m_post, m_pre, y_pre, "age_post", "education_post"), response = y_post), data = dfY)

  a_hat   <- unname(coef(fit_m)["dose"])
  b_hat   <- unname(coef(fit_y)[m_post])
  dir_hat <- unname(coef(fit_y)["dose"])
  ind_hat <- as.numeric(a_hat * b_hat)

  # standardized fits
  std_cols_m <- c(m_post, "dose", m_pre, "age_post", "education_post")
  std_cols_y <- c(y_post, "dose", m_post, m_pre, y_pre, "age_post", "education_post")

  dfM_std <- zscore_df(dfM, std_cols_m)
  dfY_std <- zscore_df(dfY, std_cols_y)

  fit_m_std <- lm(reformulate(c("dose", m_pre, "age_post", "education_post"), response = m_post), data = dfM_std)
  fit_y_std <- lm(reformulate(c("dose", m_post, m_pre, y_pre, "age_post", "education_post"), response = y_post), data = dfY_std)

  a_hat_std   <- unname(coef(fit_m_std)["dose"])
  b_hat_std   <- unname(coef(fit_y_std)[m_post])
  dir_hat_std <- unname(coef(fit_y_std)["dose"])
  ind_hat_std <- as.numeric(a_hat_std * b_hat_std)

  ind_vec <- numeric(R)
  dir_vec <- numeric(R)
  ind_vec_std <- numeric(R)
  dir_vec_std <- numeric(R)

  for (i in seq_len(R)) {
    idxM <- sample.int(nM, nM, replace = TRUE)
    idxY <- sample.int(nY, nY, replace = TRUE)

    dfiM <- dfM[idxM, , drop = FALSE]
    dfiY <- dfY[idxY, , drop = FALSE]

    fit_mi <- try(lm(reformulate(c("dose", m_pre, "age_post", "education_post"), response = m_post), data = dfiM), silent = TRUE)
    fit_yi <- try(lm(reformulate(c("dose", m_post, m_pre, y_pre, "age_post", "education_post"), response = y_post), data = dfiY), silent = TRUE)

    dfiM_std <- zscore_df(dfiM, std_cols_m)
    dfiY_std <- zscore_df(dfiY, std_cols_y)

    fit_mi_std <- try(lm(reformulate(c("dose", m_pre, "age_post", "education_post"), response = m_post), data = dfiM_std), silent = TRUE)
    fit_yi_std <- try(lm(reformulate(c("dose", m_post, m_pre, y_pre, "age_post", "education_post"), response = y_post), data = dfiY_std), silent = TRUE)

    if (inherits(fit_mi, "try-error") || inherits(fit_yi, "try-error")) {
      ind_vec[i] <- NA_real_
      dir_vec[i] <- NA_real_
    } else {
      ai <- unname(coef(fit_mi)["dose"])
      bi <- unname(coef(fit_yi)[m_post])
      ci <- unname(coef(fit_yi)["dose"])
      ind_vec[i] <- as.numeric(ai * bi)
      dir_vec[i] <- as.numeric(ci)
    }

    if (inherits(fit_mi_std, "try-error") || inherits(fit_yi_std, "try-error")) {
      ind_vec_std[i] <- NA_real_
      dir_vec_std[i] <- NA_real_
    } else {
      ai_s <- unname(coef(fit_mi_std)["dose"])
      bi_s <- unname(coef(fit_yi_std)[m_post])
      ci_s <- unname(coef(fit_yi_std)["dose"])
      ind_vec_std[i] <- as.numeric(ai_s * bi_s)
      dir_vec_std[i] <- as.numeric(ci_s)
    }
  }

  tibble(
    n_M = nM,
    n_Y = nY,
    ind_est = ind_hat,
    ind_lwr = quantile(ind_vec, 0.025, na.rm = TRUE),
    ind_upr = quantile(ind_vec, 0.975, na.rm = TRUE),
    dir_est = as.numeric(dir_hat),
    dir_lwr = quantile(dir_vec, 0.025, na.rm = TRUE),
    dir_upr = quantile(dir_vec, 0.975, na.rm = TRUE),
    n_boot  = sum(is.finite(ind_vec)),
    ind_est_std = ind_hat_std,
    ind_lwr_std = quantile(ind_vec_std, 0.025, na.rm = TRUE),
    ind_upr_std = quantile(ind_vec_std, 0.975, na.rm = TRUE),
    dir_est_std = as.numeric(dir_hat_std),
    dir_lwr_std = quantile(dir_vec_std, 0.025, na.rm = TRUE),
    dir_upr_std = quantile(dir_vec_std, 0.975, na.rm = TRUE),
    n_boot_std  = sum(is.finite(ind_vec_std))
  )
}

# ---------------- single job runner ----------------
run_one <- function(y_post, y_pre, m_post, m_pre, mediator_metric, mediator_band, mediator_region) {
  need_cols_M <- c("ID", "dose", "age_post", "education_post", m_post, m_pre)
  dM <- dat0 %>% select(any_of(need_cols_M)) %>% drop_na(all_of(need_cols_M))
  if (nrow(dM) < 30) return(NULL)

  need_cols_Y <- c("ID", "dose", "age_post", "education_post", y_post, y_pre, m_post, m_pre)
  dY <- dat0 %>% select(any_of(need_cols_Y)) %>% drop_na(all_of(need_cols_Y))
  if (nrow(dY) < 30) return(NULL)

  nM <- nrow(dM)
  nY <- nrow(dY)

  fit_m <- lm(reformulate(c("dose", m_pre, "age_post", "education_post"), response = m_post), data = dM)
  fit_y <- lm(reformulate(c("dose", m_post, m_pre, y_pre, "age_post", "education_post"), response = y_post), data = dY)

  m_info <- tidy_augment(fit_m, conf_level = 0.95)
  y_info <- tidy_augment(fit_y, conf_level = 0.95)
  sm_m <- m_info$tidy
  sm_y <- y_info$tidy
  gl_m <- m_info$glance
  gl_y <- y_info$glance

  std_cols_m <- c(m_post, "dose", m_pre, "age_post", "education_post")
  std_cols_y <- c(y_post, "dose", m_post, m_pre, y_pre, "age_post", "education_post")

  dM_std <- zscore_df(dM, std_cols_m)
  dY_std <- zscore_df(dY, std_cols_y)

  fit_m_std <- lm(reformulate(c("dose", m_pre, "age_post", "education_post"), response = m_post), data = dM_std)
  fit_y_std <- lm(reformulate(c("dose", m_post, m_pre, y_pre, "age_post", "education_post"), response = y_post), data = dY_std)

  sm_m_std <- broom::tidy(fit_m_std)
  sm_y_std <- broom::tidy(fit_y_std)

  base_cols <- tibble(
    y_post = y_post,
    y_pre  = y_pre,
    mediator_metric = mediator_metric,
    mediator_band   = mediator_band,
    mediator_region = mediator_region,
    mediator_post = m_post,
    mediator_pre  = m_pre,
    n_M = nM,
    n_Y = nY,
    n = nY
  )

  get_term_safe <- function(df, term_name) {
    out <- df %>% filter(term == !!term_name) %>% slice_head(n = 1)
    if (nrow(out) == 0) {
      return(tibble(
        estimate = NA_real_, std.error = NA_real_, statistic = NA_real_,
        p.value = NA_real_, conf.low = NA_real_, conf.high = NA_real_
      ))
    }
    out
  }

  a_row <- sm_m %>% filter(term == "dose") %>% slice_head(n = 1)
  b_row <- sm_y %>% filter(term == m_post)  %>% slice_head(n = 1)
  c_row <- sm_y %>% filter(term == "dose")  %>% slice_head(n = 1)

  a_row_std <- sm_m_std %>% filter(term == "dose") %>% slice_head(n = 1)
  b_row_std <- sm_y_std %>% filter(term == m_post)  %>% slice_head(n = 1)
  c_row_std <- sm_y_std %>% filter(term == "dose")  %>% slice_head(n = 1)

  m_Mpre <- get_term_safe(sm_m, m_pre)
  m_age  <- get_term_safe(sm_m, "age_post")
  m_edu  <- get_term_safe(sm_m, "education_post")

  y_Mpre <- get_term_safe(sm_y, m_pre)
  y_Ypre <- get_term_safe(sm_y, y_pre)
  y_age  <- get_term_safe(sm_y, "age_post")
  y_edu  <- get_term_safe(sm_y, "education_post")

  get_beta <- function(tidy_std, term) {
    x <- tidy_std %>% filter(term == !!term) %>% slice_head(n = 1)
    if (nrow(x) == 0) return(NA_real_)
    as.numeric(x$estimate)
  }

  m_Mpre_beta <- get_beta(sm_m_std, m_pre)
  m_age_beta  <- get_beta(sm_m_std, "age_post")
  m_edu_beta  <- get_beta(sm_m_std, "education_post")
  y_Mpre_beta <- get_beta(sm_y_std, m_pre)
  y_Ypre_beta <- get_beta(sm_y_std, y_pre)
  y_age_beta  <- get_beta(sm_y_std, "age_post")
  y_edu_beta  <- get_beta(sm_y_std, "education_post")

  bs <- boot_indirect_split(
    dfM = dM, dfY = dY,
    y_post = y_post, y_pre = y_pre,
    m_post = m_post, m_pre = m_pre,
    R = 2000, seed = 2025
  )

  wide_row <- tibble(
    y_post = y_post,
    y_pre  = y_pre,
    mediator_metric = mediator_metric,
    mediator_band   = mediator_band,
    mediator_region = mediator_region,
    mediator_post = m_post,
    mediator_pre  = m_pre,
    n_M = nM,
    n_Y = nY,
    n = nY,

    a_est = a_row$estimate %||% NA_real_,
    a_se  = a_row$std.error %||% NA_real_,
    a_t   = a_row$statistic %||% NA_real_,
    a_p   = a_row$p.value %||% NA_real_,
    a_ci_l= a_row$conf.low %||% NA_real_,
    a_ci_u= a_row$conf.high %||% NA_real_,
    a_est_std = a_row_std$estimate %||% NA_real_,

    b_est = b_row$estimate %||% NA_real_,
    b_se  = b_row$std.error %||% NA_real_,
    b_t   = b_row$statistic %||% NA_real_,
    b_p   = b_row$p.value %||% NA_real_,
    b_ci_l= b_row$conf.low %||% NA_real_,
    b_ci_u= b_row$conf.high %||% NA_real_,
    b_est_std = b_row_std$estimate %||% NA_real_,

    cprime_est = c_row$estimate %||% NA_real_,
    cprime_se  = c_row$std.error %||% NA_real_,
    cprime_t   = c_row$statistic %||% NA_real_,
    cprime_p   = c_row$p.value %||% NA_real_,
    cprime_ci_l= c_row$conf.low %||% NA_real_,
    cprime_ci_u= c_row$conf.high %||% NA_real_,
    cprime_est_std = c_row_std$estimate %||% NA_real_,

    m_Mpre_est = m_Mpre$estimate %||% NA_real_,
    m_Mpre_se  = m_Mpre$std.error %||% NA_real_,
    m_Mpre_t   = m_Mpre$statistic %||% NA_real_,
    m_Mpre_p   = m_Mpre$p.value %||% NA_real_,
    m_Mpre_ci_l= m_Mpre$conf.low %||% NA_real_,
    m_Mpre_ci_u= m_Mpre$conf.high %||% NA_real_,

    m_age_est = m_age$estimate %||% NA_real_,
    m_age_se  = m_age$std.error %||% NA_real_,
    m_age_t   = m_age$statistic %||% NA_real_,
    m_age_p   = m_age$p.value %||% NA_real_,
    m_age_ci_l= m_age$conf.low %||% NA_real_,
    m_age_ci_u= m_age$conf.high %||% NA_real_,

    m_edu_est = m_edu$estimate %||% NA_real_,
    m_edu_se  = m_edu$std.error %||% NA_real_,
    m_edu_t   = m_edu$statistic %||% NA_real_,
    m_edu_p   = m_edu$p.value %||% NA_real_,
    m_edu_ci_l= m_edu$conf.low %||% NA_real_,
    m_edu_ci_u= m_edu$conf.high %||% NA_real_,

    y_Mpre_est = y_Mpre$estimate %||% NA_real_,
    y_Mpre_se  = y_Mpre$std.error %||% NA_real_,
    y_Mpre_t   = y_Mpre$statistic %||% NA_real_,
    y_Mpre_p   = y_Mpre$p.value %||% NA_real_,
    y_Mpre_ci_l= y_Mpre$conf.low %||% NA_real_,
    y_Mpre_ci_u= y_Mpre$conf.high %||% NA_real_,

    y_Ypre_est = y_Ypre$estimate %||% NA_real_,
    y_Ypre_se  = y_Ypre$std.error %||% NA_real_,
    y_Ypre_t   = y_Ypre$statistic %||% NA_real_,
    y_Ypre_p   = y_Ypre$p.value %||% NA_real_,
    y_Ypre_ci_l= y_Ypre$conf.low %||% NA_real_,
    y_Ypre_ci_u= y_Ypre$conf.high %||% NA_real_,

    y_age_est = y_age$estimate %||% NA_real_,
    y_age_se  = y_age$std.error %||% NA_real_,
    y_age_t   = y_age$statistic %||% NA_real_,
    y_age_p   = y_age$p.value %||% NA_real_,
    y_age_ci_l= y_age$conf.low %||% NA_real_,
    y_age_ci_u= y_age$conf.high %||% NA_real_,

    y_edu_est = y_edu$estimate %||% NA_real_,
    y_edu_se  = y_edu$std.error %||% NA_real_,
    y_edu_t   = y_edu$statistic %||% NA_real_,
    y_edu_p   = y_edu$p.value %||% NA_real_,
    y_edu_ci_l= y_edu$conf.low %||% NA_real_,
    y_edu_ci_u= y_edu$conf.high %||% NA_real_,

    m_r2   = gl_m$r.squared %||% NA_real_,
    m_r2_adj = gl_m$adj.r.squared %||% NA_real_,
    m_sigma  = gl_m$sigma %||% NA_real_,
    m_AIC    = gl_m$AIC %||% NA_real_,
    m_BIC    = gl_m$BIC %||% NA_real_,
    m_logLik = gl_m$logLik %||% NA_real_,
    m_F      = gl_m$statistic %||% NA_real_,
    m_df1    = gl_m$df %||% NA_real_,
    m_df2    = gl_m$df.residual %||% NA_real_,
    m_F_p    = gl_m$p.value %||% NA_real_,

    y_r2   = gl_y$r.squared %||% NA_real_,
    y_r2_adj = gl_y$adj.r.squared %||% NA_real_,
    y_sigma  = gl_y$sigma %||% NA_real_,
    y_AIC    = gl_y$AIC %||% NA_real_,
    y_BIC    = gl_y$BIC %||% NA_real_,
    y_logLik = gl_y$logLik %||% NA_real_,
    y_F      = gl_y$statistic %||% NA_real_,
    y_df1    = gl_y$df %||% NA_real_,
    y_df2    = gl_y$df.residual %||% NA_real_,
    y_F_p    = gl_y$p.value %||% NA_real_,

    ind_est = bs$ind_est,
    ind_lwr = bs$ind_lwr,
    ind_upr = bs$ind_upr,
    dir_est = bs$dir_est,
    dir_lwr = bs$dir_lwr,
    dir_upr = bs$dir_upr,
    n_boot  = bs$n_boot,

    ind_est_std = bs$ind_est_std,
    ind_lwr_std = bs$ind_lwr_std,
    ind_upr_std = bs$ind_upr_std,
    dir_est_std = bs$dir_est_std,
    dir_lwr_std = bs$dir_lwr_std,
    dir_upr_std = bs$dir_upr_std,
    n_boot_std  = bs$n_boot_std,

    m_Mpre_beta = m_Mpre_beta,
    m_age_beta  = m_age_beta,
    m_edu_beta  = m_edu_beta,
    y_Mpre_beta = y_Mpre_beta,
    y_Ypre_beta = y_Ypre_beta,
    y_age_beta  = y_age_beta,
    y_edu_beta  = y_edu_beta
  )

  long_m <- sm_m %>%
    mutate(model = "M", outcome = m_post) %>%
    bind_cols(base_cols[rep(1, nrow(.)), ]) %>%
    relocate(y_post, y_pre, mediator_metric, mediator_band, mediator_region,
             mediator_post, mediator_pre, n_M, n_Y, n, model, outcome,
             term, estimate, std.error, statistic, p.value, conf.low, conf.high)

  long_y <- sm_y %>%
    mutate(model = "Y", outcome = y_post) %>%
    bind_cols(base_cols[rep(1, nrow(.)), ]) %>%
    relocate(y_post, y_pre, mediator_metric, mediator_band, mediator_region,
             mediator_post, mediator_pre, n_M, n_Y, n, model, outcome,
             term, estimate, std.error, statistic, p.value, conf.low, conf.high)

  list(long = bind_rows(long_m, long_y), wide = wide_row)
}

# ---------------- build jobs ----------------
jobs <- tidyr::expand_grid(
  y_post = y_map$y_post,
  mediators_def
) %>%
  left_join(y_map, by = "y_post") %>%
  transmute(
    y_post,
    y_pre,
    m_post = m_post,
    m_pre  = m_pre,
    mediator_metric,
    mediator_band,
    mediator_region
  )

# ---------------- parallel run ----------------
plan(multisession, workers = 10)
handlers(global = TRUE)

with_progress({
  p <- progressor(steps = nrow(jobs))
  res_list <- future_pmap(
    jobs,
    function(y_post, y_pre, m_post, m_pre, mediator_metric, mediator_band, mediator_region) {
      p()
      run_one(y_post, y_pre, m_post, m_pre, mediator_metric, mediator_band, mediator_region)
    }
  )

  res_long <- bind_rows(map(res_list, "long"))
  res_wide <- bind_rows(map(res_list, "wide"))

  readr::write_csv(res_long, file.path(out_dir, paste0(file_prefix, "mediation_terms_long.csv")))
  readr::write_csv(res_wide, file.path(out_dir, paste0(file_prefix, "mediation_full_paths_bootstrap.csv")))

  readr::write_csv(
    res_wide %>%
      select(
        y_post, y_pre, mediator_metric, mediator_band, mediator_region,
        mediator_post, mediator_pre, n, n_M, n_Y,
        a_est, a_p, b_est, b_p, cprime_est, cprime_p,
        ind_est, ind_lwr, ind_upr, dir_est, dir_lwr, dir_upr, n_boot,
        a_est_std, b_est_std, cprime_est_std,
        ind_est_std, ind_lwr_std, ind_upr_std, dir_est_std, dir_lwr_std, dir_upr_std
      ),
    file.path(out_dir, paste0(file_prefix, "mediation_dose_Mpost_Ypost_bootstrap.csv"))
  )

  sig_indirect_ci <- res_wide %>%
    filter(is.finite(ind_lwr), is.finite(ind_upr)) %>%
    mutate(ind_sig95 = (ind_lwr > 0 | ind_upr < 0)) %>%
    filter(ind_sig95) %>%
    arrange(y_post, mediator_metric, mediator_band, mediator_region)
  readr::write_csv(sig_indirect_ci, file.path(out_dir, paste0(file_prefix, "mediation_indirect_significant_95CI.csv")))

  sig_ab <- res_wide %>%
    filter(!is.na(a_p), !is.na(b_p)) %>%
    filter(a_p < 0.05, b_p < 0.05) %>%
    arrange(y_post, mediator_metric, mediator_band, mediator_region)
  readr::write_csv(sig_ab, file.path(out_dir, paste0(file_prefix, "mediation_indirect_sig_by_ab_paths.csv")))

  recommended <- bind_rows(
    sig_indirect_ci %>% mutate(flag = "CI_no_zero"),
    sig_ab         %>% mutate(flag = "a&b_p<0.05")
  ) %>%
    distinct(y_post, mediator_metric, mediator_band, mediator_region, .keep_all = TRUE) %>%
    arrange(y_post, mediator_metric, mediator_band, mediator_region)
  readr::write_csv(recommended, file.path(out_dir, paste0(file_prefix, "mediation_recommended_significant_list.csv")))

  keys <- c("y_post", "y_pre", "mediator_metric", "mediator_band", "mediator_region", "mediator_post", "mediator_pre")
  keep_cols <- c(keys, "n", "n_M", "n_Y", "a_est", "a_p", "b_est", "b_p", "cprime_est", "cprime_p",
                 "ind_est", "ind_lwr", "ind_upr", "dir_est", "dir_lwr", "dir_upr", "n_boot")

  ci_s <- sig_indirect_ci %>% select(any_of(keep_cols)) %>% mutate(ind_sig95 = TRUE)
  ab_s <- sig_ab %>% select(any_of(keep_cols)) %>% mutate(ab_sig = TRUE)
  merged <- full_join(ci_s, ab_s, by = keys, suffix = c(".ci", ".ab"))
  coalesce2 <- function(x_ci, x_ab) dplyr::coalesce(x_ci, x_ab)

  recommended_both <- merged %>%
    mutate(
      ind_sig95 = if_else(is.na(ind_sig95), FALSE, ind_sig95),
      ab_sig    = if_else(is.na(ab_sig), FALSE, ab_sig),
      n          = coalesce2(n.ci, n.ab),
      n_M        = coalesce2(n_M.ci, n_M.ab),
      n_Y        = coalesce2(n_Y.ci, n_Y.ab),
      a_est      = coalesce2(a_est.ci, a_est.ab),
      a_p        = coalesce2(a_p.ci, a_p.ab),
      b_est      = coalesce2(b_est.ci, b_est.ab),
      b_p        = coalesce2(b_p.ci, b_p.ab),
      cprime_est = coalesce2(cprime_est.ci, cprime_est.ab),
      cprime_p   = coalesce2(cprime_p.ci, cprime_p.ab),
      ind_est    = coalesce2(ind_est.ci, ind_est.ab),
      ind_lwr    = coalesce2(ind_lwr.ci, ind_lwr.ab),
      ind_upr    = coalesce2(ind_upr.ci, ind_upr.ab),
      dir_est    = coalesce2(dir_est.ci, dir_est.ab),
      dir_lwr    = coalesce2(dir_lwr.ci, dir_lwr.ab),
      dir_upr    = coalesce2(dir_upr.ci, dir_upr.ab),
      n_boot     = coalesce2(n_boot.ci, n_boot.ab)
    ) %>%
    mutate(
      flag = case_when(
        ind_sig95 & ab_sig ~ "both",
        ind_sig95          ~ "CI_no_zero",
        ab_sig             ~ "a&b_p<0.05",
        TRUE               ~ NA_character_
      )
    ) %>%
    filter(ind_sig95 | ab_sig) %>%
    select(all_of(keys),
           n, n_M, n_Y, a_est, a_p, b_est, b_p, cprime_est, cprime_p,
           ind_est, ind_lwr, ind_upr, dir_est, dir_lwr, dir_upr, n_boot,
           ind_sig95, ab_sig, flag) %>%
    arrange(y_post, mediator_metric, mediator_band, mediator_region)
  readr::write_csv(recommended_both, file.path(out_dir, paste0(file_prefix, "mediation_recommended_significant_list_with_both.csv")))

  widen_terms <- function(df, model_tag) {
    df %>%
      mutate(term = paste0(model_tag, "__", term)) %>%
      select(y_post, y_pre, mediator_metric, mediator_band, mediator_region,
             mediator_post, mediator_pre, n,
             term, estimate, std.error, statistic, p.value, conf.low, conf.high) %>%
      pivot_wider(
        names_from = term,
        values_from = c(estimate, std.error, statistic, p.value, conf.low, conf.high),
        names_sep = "."
      )
  }

  build_all_paths_unstd <- function(keys_tbl) {
    long_sig <- res_long %>%
      inner_join(keys_tbl,
                 by = c("y_post", "y_pre", "mediator_metric", "mediator_band", "mediator_region",
                        "mediator_post", "mediator_pre", "n"))
    wide_M <- long_sig %>% filter(model == "M") %>% widen_terms("M")
    wide_Y <- long_sig %>% filter(model == "Y") %>% widen_terms("Y")
    list(wide_M = wide_M, wide_Y = wide_Y)
  }

  get_std_tidy_for_job <- function(job_row) {
    y_post <- job_row$y_post
    y_pre  <- job_row$y_pre
    m_post <- job_row$mediator_post
    m_pre  <- job_row$mediator_pre

    need_cols_M <- c("ID", "dose", "age_post", "education_post", m_post, m_pre)
    dM <- dat0 %>% select(any_of(need_cols_M)) %>% drop_na(all_of(need_cols_M))
    if (nrow(dM) < 30) return(NULL)

    need_cols_Y <- c("ID", "dose", "age_post", "education_post", y_post, y_pre, m_post, m_pre)
    dY <- dat0 %>% select(any_of(need_cols_Y)) %>% drop_na(all_of(need_cols_Y))
    if (nrow(dY) < 30) return(NULL)

    std_cols_m <- c(m_post, "dose", m_pre, "age_post", "education_post")
    std_cols_y <- c(y_post, "dose", m_post, m_pre, y_pre, "age_post", "education_post")

    dM_std <- zscore_df(dM, std_cols_m)
    dY_std <- zscore_df(dY, std_cols_y)

    fit_m_std <- lm(reformulate(c("dose", m_pre, "age_post", "education_post"), response = m_post), data = dM_std)
    fit_y_std <- lm(reformulate(c("dose", m_post, m_pre, y_pre, "age_post", "education_post"), response = y_post), data = dY_std)

    sm_m_std <- broom::tidy(fit_m_std, conf.int = TRUE)
    sm_y_std <- broom::tidy(fit_y_std, conf.int = TRUE)

    base_cols <- tibble(
      y_post = y_post,
      y_pre  = y_pre,
      mediator_post = m_post,
      mediator_pre  = m_pre,
      mediator_metric = job_row$mediator_metric,
      mediator_band   = job_row$mediator_band,
      mediator_region = job_row$mediator_region,
      n = nrow(dY)
    )

    long_m_std <- sm_m_std %>% mutate(model = "M") %>% bind_cols(base_cols[rep(1, nrow(.)), ])
    long_y_std <- sm_y_std %>% mutate(model = "Y") %>% bind_cols(base_cols[rep(1, nrow(.)), ])
    bind_rows(long_m_std, long_y_std)
  }

  build_all_paths_std <- function(keys_tbl) {
    sig_jobs <- keys_tbl %>%
      select(y_post, y_pre, mediator_metric, mediator_band, mediator_region, mediator_post, mediator_pre, n) %>%
      distinct()

    std_tidy_list <- purrr::pmap(sig_jobs, function(y_post, y_pre, mediator_metric, mediator_band, mediator_region, mediator_post, mediator_pre, n) {
      get_std_tidy_for_job(tibble(
        y_post = y_post, y_pre = y_pre,
        mediator_metric = mediator_metric,
        mediator_band = mediator_band,
        mediator_region = mediator_region,
        mediator_post = mediator_post,
        mediator_pre = mediator_pre
      ))
    }) %>% purrr::compact()

    if (length(std_tidy_list) == 0) return(list(wide_M = NULL, wide_Y = NULL))
    long_sig_std <- bind_rows(std_tidy_list)
    wide_M_std <- long_sig_std %>% filter(model == "M") %>% widen_terms("M")
    wide_Y_std <- long_sig_std %>% filter(model == "Y") %>% widen_terms("Y")
    list(wide_M = wide_M_std, wide_Y = wide_Y_std)
  }

  keep_main_unstd <- c(
    "y_post", "y_pre", "mediator_metric", "mediator_band", "mediator_region",
    "mediator_post", "mediator_pre", "n", "n_M", "n_Y",
    "a_est", "a_se", "a_t", "a_p", "a_ci_l", "a_ci_u",
    "b_est", "b_se", "b_t", "b_p", "b_ci_l", "b_ci_u",
    "cprime_est", "cprime_se", "cprime_t", "cprime_p", "cprime_ci_l", "cprime_ci_u",
    "ind_est", "ind_lwr", "ind_upr", "dir_est", "dir_lwr", "dir_upr", "n_boot",
    "m_Mpre_beta", "m_age_beta", "m_edu_beta",
    "y_Mpre_beta", "y_Ypre_beta", "y_age_beta", "y_edu_beta"
  )

  keep_main_std <- c(
    "y_post", "y_pre", "mediator_metric", "mediator_band", "mediator_region",
    "mediator_post", "mediator_pre", "n",
    "a_est_std", "b_est_std", "cprime_est_std",
    "ind_est_std", "ind_lwr_std", "ind_upr_std", "dir_est_std", "dir_lwr_std", "dir_upr_std", "n_boot_std"
  )

  assemble_full_table <- function(keys_tbl, which = c("unstd", "std"), base_filename) {
    which <- match.arg(which)
    if (nrow(keys_tbl) == 0) return(invisible(NULL))

    if (which == "unstd") {
      ap <- build_all_paths_unstd(keys_tbl)
      main_block <- res_wide %>%
        inner_join(keys_tbl,
                   by = c("y_post", "y_pre", "mediator_metric", "mediator_band", "mediator_region",
                          "mediator_post", "mediator_pre", "n")) %>%
        select(any_of(keep_main_unstd))

      full_tab <- main_block %>%
        left_join(ap$wide_M, by = c("y_post", "y_pre", "mediator_metric", "mediator_band", "mediator_region",
                                    "mediator_post", "mediator_pre", "n")) %>%
        left_join(ap$wide_Y, by = c("y_post", "y_pre", "mediator_metric", "mediator_band", "mediator_region",
                                    "mediator_post", "mediator_pre", "n"))

      readr::write_csv(full_tab, file.path(out_dir, paste0(file_prefix, base_filename, "_all_terms_unstd.csv")))
    } else {
      ap <- build_all_paths_std(keys_tbl)
      main_block <- res_wide %>%
        inner_join(keys_tbl,
                   by = c("y_post", "y_pre", "mediator_metric", "mediator_band", "mediator_region",
                          "mediator_post", "mediator_pre", "n")) %>%
        select(any_of(keep_main_std))

      full_tab <- main_block
      if (!is.null(ap$wide_M) && !is.null(ap$wide_Y)) {
        full_tab <- full_tab %>%
          left_join(ap$wide_M, by = c("y_post", "y_pre", "mediator_metric", "mediator_band", "mediator_region",
                                      "mediator_post", "mediator_pre", "n")) %>%
          left_join(ap$wide_Y, by = c("y_post", "y_pre", "mediator_metric", "mediator_band", "mediator_region",
                                      "mediator_post", "mediator_pre", "n"))
      }
      readr::write_csv(full_tab, file.path(out_dir, paste0(file_prefix, base_filename, "_all_terms_std.csv")))
    }
  }

  sig_ci_keys <- res_wide %>%
    filter(is.finite(ind_lwr), is.finite(ind_upr)) %>%
    mutate(ind_sig95 = (ind_lwr > 0 | ind_upr < 0)) %>%
    filter(ind_sig95) %>%
    select(y_post, y_pre, mediator_metric, mediator_band, mediator_region,
           mediator_post, mediator_pre, n) %>%
    distinct()
  assemble_full_table(sig_ci_keys, "unstd", "mediation_indirect_significant_95CI")
  assemble_full_table(sig_ci_keys, "std",   "mediation_indirect_significant_95CI")

  sig_ab_keys <- res_wide %>%
    filter(!is.na(a_p), !is.na(b_p)) %>%
    filter(a_p < 0.05, b_p < 0.05) %>%
    select(y_post, y_pre, mediator_metric, mediator_band, mediator_region,
           mediator_post, mediator_pre, n) %>%
    distinct()
  assemble_full_table(sig_ab_keys, "unstd", "mediation_indirect_sig_by_ab_paths")
  assemble_full_table(sig_ab_keys, "std",   "mediation_indirect_sig_by_ab_paths")

  rec_keys <- recommended_both %>%
    select(y_post, y_pre, mediator_metric, mediator_band, mediator_region,
           mediator_post, mediator_pre, n) %>%
    distinct()
  assemble_full_table(rec_keys, "unstd", "mediation_recommended_significant_list_with_both")
  assemble_full_table(rec_keys, "std",   "mediation_recommended_significant_list_with_both")

  message("Output dir: ", normalizePath(out_dir, mustWork = FALSE))

  # ================= original TOPZ8 outputs retained =================
  topz8_keys <- tribble(
    ~y_post, ~mediator_metric, ~mediator_band, ~mediator_region,
    "pronesupinelh_post",         "topz_rate",   "HighAlpha",     "Temporal",
    "rigidlowerl_post",           "topz_rate",   "HighAlpha",     "Temporal",
    "rigidupperl_post",           "topz_rate",   "HighAlpha",     "Temporal",
    "taplf_post",                 "topz_rate",   "HighAlpha",     "Temporal",
    "wmSequence_a1_1_post",       "topz_rate",   "HighAlpha",     "Temporal",
    "restingtremorlowerr_post",   "topz_rate",   "HighAlpha",     "Occipital",
    "speech_post",                "topz_rate",   "HighAlpha",     "Occipital",
    "tiempototal2_a1_1_post",     "topz_rate",   "Beta",          "Temporal"
  )

  topz8_out <- res_wide %>%
    semi_join(topz8_keys,
              by = c("y_post", "mediator_metric", "mediator_band", "mediator_region")) %>%
    arrange(mediator_band, mediator_region, y_post)

  if (nrow(topz8_out) != nrow(topz8_keys)) {
    missing8 <- topz8_keys %>%
      anti_join(res_wide,
                by = c("y_post", "mediator_metric", "mediator_band", "mediator_region"))
    warning(
      "TopZ8: Not all 8 rows were found in res_wide.\nMissing rows:\n",
      paste(capture.output(print(missing8)), collapse = "\n")
    )
  }

  readr::write_csv(topz8_out, file.path(out_dir, paste0(file_prefix, "TOPZ8_full_paths_bootstrap.csv")))

  topz8_slim <- topz8_out %>%
    select(
      y_post, y_pre,
      mediator_metric, mediator_band, mediator_region,
      mediator_post, mediator_pre,
      n, n_M, n_Y,
      a_est, a_p, a_est_std,
      b_est, b_p, b_est_std,
      cprime_est, cprime_p, cprime_est_std,
      ind_est, ind_lwr, ind_upr,
      ind_est_std, ind_lwr_std, ind_upr_std,
      dir_est, dir_lwr, dir_upr,
      dir_est_std, dir_lwr_std, dir_upr_std,
      n_boot, n_boot_std
    )

  readr::write_csv(topz8_slim, file.path(out_dir, paste0(file_prefix, "TOPZ8_key_results.csv")))

  # ================= retained diagnostic extraction =================
  cat("\n========== EXTRACTING GENUINE TopZ BOTH RESULTS ==========\n")

  topz_both_genuine <- res_wide %>%
    filter(mediator_metric == "topz_rate") %>%
    filter(!is.na(a_p), !is.na(b_p)) %>%
    filter(a_p < 0.05, b_p < 0.05) %>%
    filter(is.finite(ind_lwr_std), is.finite(ind_upr_std)) %>%
    filter(ind_lwr_std > 0 | ind_upr_std < 0) %>%
    arrange(mediator_band, mediator_region, y_post)

  cat("Found", nrow(topz_both_genuine), "TopZ mediations meeting genuine BOTH criteria\n")

  topz_summary <- topz_both_genuine %>%
    count(mediator_band, mediator_region, name = "n_outcomes") %>%
    arrange(desc(n_outcomes))

  cat("\nBand-Region distribution:\n")
  print(topz_summary, n = Inf)

  readr::write_csv(topz_both_genuine, file.path(out_dir, paste0(file_prefix, "TOPZ_BOTH_Genuine_full_paths_bootstrap.csv")))

  topz_both_key <- topz_both_genuine %>%
    select(
      y_post, y_pre,
      mediator_metric, mediator_band, mediator_region,
      mediator_post, mediator_pre,
      n, n_M, n_Y,
      a_est, a_p, a_est_std,
      b_est, b_p, b_est_std,
      cprime_est, cprime_p, cprime_est_std,
      ind_est, ind_lwr, ind_upr,
      ind_est_std, ind_lwr_std, ind_upr_std,
      dir_est, dir_lwr, dir_upr,
      dir_est_std, dir_lwr_std, dir_upr_std,
      n_boot, n_boot_std
    )
  readr::write_csv(topz_both_key, file.path(out_dir, paste0(file_prefix, "TOPZ_BOTH_Genuine_key_results.csv")))

  topz_both_motor <- topz_both_genuine %>%
    filter(!str_detect(y_post, regex("wm|stroop|tiempototal|sequence|memoria", ignore_case = TRUE))) %>%
    select(
      y_post, y_pre,
      mediator_metric, mediator_band, mediator_region,
      mediator_post, mediator_pre,
      n, a_est, a_p, a_est_std,
      b_est, b_p, b_est_std,
      cprime_est, cprime_p, cprime_est_std,
      ind_est_std, ind_lwr_std, ind_upr_std,
      n_boot_std
    )
  readr::write_csv(topz_both_motor, file.path(out_dir, paste0(file_prefix, "TOPZ_BOTH_Genuine_motor_only.csv")))

  cat("\n✓ Output files generated:\n")
  cat("  - TOPZ_BOTH_Genuine_full_paths_bootstrap.csv (all columns)\n")
  cat("  - TOPZ_BOTH_Genuine_key_results.csv (for Figure 3A)\n")
  cat("  - TOPZ_BOTH_Genuine_motor_only.csv (motor outcomes only)\n")

  cat("\n========== COMPARING WITH ORIGINAL TOPZ8 ==========\n")

  topz8_keys_original <- tribble(
    ~y_post, ~mediator_metric, ~mediator_band, ~mediator_region,
    "pronesupinelh_post",         "topz_rate",   "HighAlpha",     "Temporal",
    "rigidlowerl_post",           "topz_rate",   "HighAlpha",     "Temporal",
    "rigidupperl_post",           "topz_rate",   "HighAlpha",     "Temporal",
    "taplf_post",                 "topz_rate",   "HighAlpha",     "Temporal",
    "wmSequence_a1_1_post",       "topz_rate",   "HighAlpha",     "Temporal",
    "restingtremorlowerr_post",   "topz_rate",   "HighAlpha",     "Occipital",
    "speech_post",                "topz_rate",   "HighAlpha",     "Occipital",
    "tiempototal2_a1_1_post",     "topz_rate",   "Beta",          "Temporal"
  )

  topz8_validation <- res_wide %>%
    inner_join(topz8_keys_original,
               by = c("y_post", "mediator_metric", "mediator_band", "mediator_region")) %>%
    mutate(
      a_sig = a_p < 0.05,
      b_sig = b_p < 0.05,
      ci_no_zero_std = is.finite(ind_lwr_std) & is.finite(ind_upr_std) &
        (ind_lwr_std > 0 | ind_upr_std < 0),
      meets_both = a_sig & b_sig & ci_no_zero_std
    ) %>%
    select(y_post, mediator_band, mediator_region,
           a_p, b_p, ind_lwr_std, ind_upr_std, meets_both)

  cat("\nOriginal TOPZ8 validation:\n")
  print(topz8_validation, n = Inf)

  cat("\nSummary:\n")
  cat("  Passed BOTH:", sum(topz8_validation$meets_both), "/ 8\n")
  cat("  Failed BOTH:", sum(!topz8_validation$meets_both), "/ 8\n")

  if (any(!topz8_validation$meets_both)) {
    cat("\n⚠️ WARNING: The following rows from original TOPZ8 do NOT meet BOTH criteria:\n")
    topz8_validation %>%
      filter(!meets_both) %>%
      select(y_post, mediator_band, mediator_region) %>%
      print(n = Inf)
  }

  cat("\n========== EXTRACTION COMPLETE ==========\n")
})
