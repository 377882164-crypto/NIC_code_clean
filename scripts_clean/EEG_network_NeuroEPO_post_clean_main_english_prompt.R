# ============================ Clean EEG network script for NIC mainline ============================
# Purpose:
#   Condition-specific PCA (scree/elbow + manual n_use) -> loading > 0.05 -> EBICglasso network
#   -> spinglass community export for downstream TOPZ pipeline.
# Notes:
#   This clean mainline script intentionally keeps only the NIC-relevant network backbone.
#   Removed from mainline: network beautification plots, legends, convex hull plots,
#   centrality exports, small-world/global/local metrics, and other non-TOPZ side outputs.

packages <- c("ggplot2", "psych", "qgraph", "igraph", "Matrix", "dplyr", "tibble", "readr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

find_project_root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    has_raw <- dir.exists(file.path(cur, "01_raw_data"))
    has_results <- dir.exists(file.path(cur, "05_results"))
    if (has_raw && has_results) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("Project root not found. Please run from inside NIC_reproduction_from_raw.")
}

project_root <- find_project_root()
condition_name <- "NeuroEPO_post"
time_label <- "post"
data_file <- file.path(project_root, "01_raw_data", "wide_dose_Hz_Neuro.csv")
out_dir <- file.path(project_root, "05_results", "network", paste0(condition_name, "_output"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("Condition: ", condition_name)
message("Input: ", data_file)
message("Output folder: ", out_dir)

# ============================ Read data and select condition-specific variables ============================
data <- read.csv(data_file, stringsAsFactors = FALSE, check.names = FALSE)
cols_all <- colnames(data)
cols_selected <- cols_all[grepl("Hz", cols_all, ignore.case = TRUE) & grepl(time_label, cols_all, ignore.case = TRUE)]
if (length(cols_selected) == 0) stop("No columns matched Hz + ", time_label, ".")

selected_data <- data[, cols_selected, drop = FALSE]
selected_data_num <- as.data.frame(lapply(selected_data, function(x) suppressWarnings(as.numeric(as.character(x)))))
selected_data_num <- na.omit(selected_data_num)
if (nrow(selected_data_num) < 3) stop("Too few complete rows after NA removal.")
selected_data_scaled <- scale(selected_data_num)

# ============================ PCA and scree/elbow ============================
cat("Running PCA and generating scree plots...\n")
pca_result <- prcomp(selected_data_scaled, center = TRUE, scale. = FALSE)
eigenvalues <- pca_result$sdev^2
n_pc <- length(eigenvalues)
var_explained <- eigenvalues / sum(eigenvalues)
cum_explained <- cumsum(var_explained)
num_pc80 <- which(cum_explained >= 0.8)[1]
cat(sprintf("Reference PC count at 80%% cumulative variance: %d\n", num_pc80))

x_pc <- seq_len(n_pc)
y_pc <- eigenvalues
line_start <- c(x_pc[1], y_pc[1])
line_end <- c(x_pc[n_pc], y_pc[n_pc])
point_line_dist <- function(x0, y0, x1, y1, x2, y2) {
  abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / sqrt((y2 - y1)^2 + (x2 - x1)^2)
}
distances <- mapply(point_line_dist, x_pc, y_pc,
                    MoreArgs = list(x1 = line_start[1], y1 = line_start[2], x2 = line_end[1], y2 = line_end[2]))
elbow <- which.max(distances)
cat(sprintf("Suggested PC count by scree elbow: %d\n", elbow))

png(file.path(out_dir, "PCA_scree_plots.png"), width = 1200, height = 600, res = 120)
par(mfrow = c(1, 2))
plot(cum_explained, type = "b", xlab = "Number of principal components", ylab = "Cumulative explained variance", main = "Cumulative explained variance")
abline(h = 0.8, col = "red", lty = 2)
points(num_pc80, cum_explained[num_pc80], col = "blue", pch = 19, cex = 1.5)
text(num_pc80, cum_explained[num_pc80], labels = paste0(num_pc80), pos = 4, col = "blue")
plot(x_pc, y_pc, type = "b", pch = 19, xlab = "Number of principal components", ylab = "Eigenvalue", main = "Scree plot")
points(elbow, y_pc[elbow], col = "red", pch = 19, cex = 1.5)
text(elbow, y_pc[elbow], labels = paste0("Elbow: ", elbow), pos = 4, col = "red")
abline(v = elbow, lty = 2, col = "red")
par(mfrow = c(1, 1))
dev.off()

cat("\n=======================\n")
cat(sprintf("Reference PC count (80%% cumulative variance): %d\n", num_pc80))
cat(sprintf("Suggested PC count (scree elbow): %d\n", elbow))
cat("Choose the number of principal components to retain:\n")
cat("  1. Enter", elbow, "to use the scree elbow recommendation\n")
cat("  2. Or enter a custom number of principal components\n")

cat("=======================\n")

n_use <- as.integer(readline(prompt = sprintf("Enter the number of principal components to retain (recommended elbow = %d): ", elbow)))
if (is.na(n_use) || n_use < 1 || n_use > n_pc) {
  cat(sprintf("Invalid input. Falling back to the 80%% cumulative-variance reference count (%d).\n", num_pc80))
  n_use <- num_pc80
}
cat(sprintf("PC count used for filtering: %d\n", n_use))

# ============================ loading > 0.05 variable retention ============================
loadings_abs <- abs(pca_result$rotation[, 1:n_use, drop = FALSE])
selected_vars_pca_idx <- which(apply(loadings_abs, 1, function(x) any(x > 0.05)))
selected_names <- colnames(selected_data_scaled)[selected_vars_pca_idx]
selected_data_pca <- selected_data_scaled[, selected_names, drop = FALSE]
cat(sprintf("Number of variables retained for downstream clustering: %d\n\n", length(selected_names)))

pca_filtered_data_path <- file.path(out_dir, "PCA_filtered_data_for_tSNE.csv")
write.csv(selected_data_pca, pca_filtered_data_path, row.names = FALSE)
message("Saved PCA-filtered variable matrix to: ", pca_filtered_data_path)

# ============================ Regularized correlation network ============================
dat_raw <- read.csv(pca_filtered_data_path, check.names = FALSE, stringsAsFactors = FALSE)
dat_num <- dat_raw %>% mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))
na_rows <- sum(!complete.cases(dat_num))
if (na_rows > 0) message("Detected ", na_rows, " row(s) with NA; removed before network estimation.")
dat_num <- stats::na.omit(dat_num)
dat_scaled <- scale(dat_num)
orig_names <- colnames(dat_scaled)

parse_channel_freq <- function(x) {
  x1 <- sub("^X+", "", x)
  x1 <- sub(paste0("(?i)(?:_)?(?:", time_label, ")$"), "", x1, perl = TRUE)
  m <- regexpr("(?i)\\d+(?:\\.\\d+)?(?=\\s*Hz)", x1, perl = TRUE)
  if (m[1] > 0) {
    freq_val <- substr(x1, m[1], m[1] + attr(m, "match.length") - 1)
    channel_raw <- gsub("(?i)\\d+(?:\\.\\d+)?\\s*Hz", "", x1, perl = TRUE)
  } else {
    m2 <- regexpr("\\d+(?:\\.\\d+)?", x1, perl = TRUE)
    if (m2[1] > 0) {
      freq_val <- substr(x1, m2[1], m2[1] + attr(m2, "match.length") - 1)
      channel_raw <- paste0(substr(x1, 1, m2[1] - 1), substr(x1, m2[1] + attr(m2, "match.length"), nchar(x1)))
    } else {
      freq_val <- "UNK"
      channel_raw <- x1
    }
  }
  channel <- gsub("[^A-Za-z0-9]", "", channel_raw)
  channel <- sub("^[0-9]+", "", channel)
  if (!nzchar(channel)) channel <- "UNK"
  list(channel = channel, freq = freq_val)
}

parsed <- lapply(orig_names, parse_channel_freq)
chan_vec <- vapply(parsed, `[[`, character(1), "channel")
freq_vec <- vapply(parsed, `[[`, character(1), "freq")
pretty_names <- paste(chan_vec, freq_vec, sep = "_")
if (anyDuplicated(pretty_names)) {
  dup_index <- ave(pretty_names, pretty_names, FUN = seq_along)
  pretty_names <- ifelse(dup_index == 1, pretty_names, paste0(pretty_names, "_dup", dup_index))
}

cor_mat <- cor(dat_scaled)
cor_mat_smooth <- psych::cor.smooth(cor_mat)
eig_vals <- eigen(cor_mat_smooth, symmetric = TRUE)$values
cor_mat_pd <- if (all(eig_vals > 1e-8)) cor_mat_smooth else {
  message("Correlation matrix is not strictly positive definite; applying nearPD correction.")
  as.matrix(Matrix::nearPD(cor_mat_smooth)$mat)
}

set.seed(123)
network <- EBICglasso(cor_mat_pd, n = nrow(dat_scaled), gamma = 0.5)
dimnames(network) <- list(pretty_names, pretty_names)

nw_df <- as.data.frame(network, check.names = FALSE) %>% rownames_to_column(var = "node")
readr::write_csv(nw_df, file.path(out_dir, "network_weights_PCAfiltered.csv"))
message("Saved regularized network weight matrix to: ", file.path(out_dir, "network_weights_PCAfiltered.csv"))

# ============================ Spinglass community export (for downstream TOPZ) ============================
g <- igraph::graph_from_adjacency_matrix(network, mode = "undirected", weighted = TRUE, diag = FALSE)
if (!igraph::is.connected(g)) {
  comps <- igraph::components(g)
  largest_id <- which.max(comps$csize)
  v_largest <- igraph::V(g)[comps$membership == largest_id]
  g_largest <- igraph::induced_subgraph(g, v_largest)
  set.seed(123)
  comm_largest <- igraph::cluster_spinglass(g_largest, weights = igraph::E(g_largest)$weight, spins = 25, implementation = "neg")
  membership_vec <- rep(NA_integer_, igraph::vcount(g))
  membership_vec[comps$membership == largest_id] <- membership(comm_largest)
  names(membership_vec) <- igraph::V(g)$name
} else {
  set.seed(123)
  comm <- igraph::cluster_spinglass(g, weights = igraph::E(g)$weight, spins = 25, implementation = "neg")
  membership_vec <- membership(comm)
}

breaks_fg <- c(0.5, 4, 8, 10, 13, 30)
labels_fg <- c("Delta 0.5--4", "Theta 4--8", "Low-Alpha 8--10", "High-Alpha 10--13", "Beta 13--30")
nodes_export <- colnames(network)
comm_aligned <- as.integer(membership_vec[nodes_export])
Channel_export <- chan_vec[match(nodes_export, pretty_names)]
Freq_raw_export <- freq_vec[match(nodes_export, pretty_names)]
FreqGroup_export <- as.character(cut(suppressWarnings(as.numeric(Freq_raw_export)),
                                     breaks = breaks_fg,
                                     labels = labels_fg,
                                     include.lowest = TRUE,
                                     right = FALSE))

detail_df <- data.frame(
  CommunityID = ifelse(is.na(comm_aligned), NA_character_, paste0("C", comm_aligned)),
  CommunityIndex = comm_aligned,
  Node = nodes_export,
  Channel = Channel_export,
  Frequency = suppressWarnings(as.numeric(Freq_raw_export)),
  FreqGroup = FreqGroup_export,
  stringsAsFactors = FALSE
)
detail_df <- detail_df[order(detail_df$CommunityIndex, detail_df$Node), ]
summary_df <- aggregate(Node ~ CommunityID + CommunityIndex, data = detail_df[!is.na(detail_df$CommunityIndex), ], FUN = length)
if (nrow(summary_df) > 0) names(summary_df)[names(summary_df) == "Node"] <- "N_Variables"
summary_df <- summary_df[order(summary_df$CommunityIndex), ]

readr::write_csv(detail_df, file.path(out_dir, "Community_Variables.csv"))
readr::write_csv(summary_df, file.path(out_dir, "Community_Summary.csv"))
message("Exported community detail and summary tables to: ", out_dir)

message("Completed clean NIC mainline EEG network pipeline: ", condition_name)
