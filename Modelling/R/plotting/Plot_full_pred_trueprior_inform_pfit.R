#!/usr/bin/env Rscript

# ============================================================
# Plot predictions from stitched RDS using an existing mapping0
# - mapping file: data/species_map.csv  (columns: species_id,species)
# - drops truth-only species
# - aligns by per-species nonzero counts (STRICT/TRUNCATE)
# - sweep_aggregation: "all" | "mean" | "mode" | "map_normconst" | "wmean_normconst"
# - OUTPUTS (7 files + extra overlays/axes):
#     1) Errors_by_sample_all.jpeg
#     2) Errors_by_sample_le20.jpeg
#     3) Predict_fullset_slab_scatter_<agg>.jpeg
#     4) Predict_fullset_summary_<agg>.jpeg                         (overlay)
#     5) Predict_fullset_summary_cutoff_<agg>.jpeg                  (overlay)
#     6) Errors_by_species_true_ge5_<agg>.jpeg
#     7) Errors_by_species_true_le5_<agg>.jpeg
#     8) Predict_fullset_summary_LOGXY_<agg>.jpeg                   (overlay, log–log)
#     9) Predict_fullset_summary_cutoff_LOGXY_<agg>.jpeg            (overlay, log–log)
#    10) Predict_fullset_summary_LOGX_<agg>.jpeg                    (overlay, log-x only)
#    11) Predict_fullset_summary_cutoff_LOGX_<agg>.jpeg             (overlay, log-x only)
# ============================================================

# -------- CONFIG (edit paths if needed) --------
library(here)
repo_root <- here::here()

rds_path        <- file.path(repo_root, "Results/step2_predictions/stitched/fullset/inform_pfit/stiched_pred_fullset.rds")
counts_path     <- file.path(repo_root, "Data/cleaned_nochimera_MATCHED_cluster_counts_ELA001_HOMOGEN.csv")
gt_path         <- file.path(repo_root, "Data/Ground_truth_adjustedorder.tsv")
species_map_csv <- file.path(repo_root, "Data/species_map.csv")

# Alignment: "strict" or "truncate"
alignment_mode     <- "strict"

# Sweep aggregation (default "all" to match prior visuals)
sweep_aggregation  <- "wmode_normconst"   # "all" | "mean" | "mode" | "map_normconst" | "wmean_normconst" | "wmode_normconst"

out_dir          <- file.path(repo_root, "Results/plots/fullset/inform_pfit")

# -------- LIBRARIES --------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(scales)  # percent axis labels + pseudo_log_trans
})

# -------- PRIOR (model-faithful: one p per species, many n per species) --------
sample_p <- function(n_spp, seed = 1) {
  set.seed(seed)
  u <- runif(n_spp)
  p <- ifelse(u < 0.65, 0.95, rnorm(n_spp, 0.5, 0.1))  # 65% spike at p=0.95; rest ~N(0.5, 0.1)
  p[p < 0.001] <- 0.001; p[p > 0.999] <- 0.999
  p
}

# Species-locked prior predictive (optional cap for plotting tails)
prior_counts_by_species <- function(df_rows, cap = 200L, seed = 42) {
  set.seed(seed)
  spp_sizes <- df_rows %>% count(species, name = "n_rows")
  ps <- sample_p(nrow(spp_sizes), seed = seed)
  spp_draws <- tibble::tibble(species = spp_sizes$species, p = ps)
  draws <- df_rows %>%
    select(species) %>%
    left_join(spp_draws, by = "species") %>%
    mutate(n = stats::rgeom(n(), prob = p) + 1L) %>%
    mutate(n = if (!is.null(cap)) pmin(n, cap) else n) %>%
    pull(n)
  draws
}

# -------- 1) Load stitched predictions & build long table --------
res <- readRDS(rds_path)
stopifnot(is.list(res), all(c("n","offsets","lengths","species_indices") %in% names(res)))

get_sweep_n <- function(x, i) {
  off <- x$offsets[i]; len <- x$lengths[i]
  if (is.na(off) || is.na(len) || len <= 0) return(numeric(0))
  x$n[(off + 1):(off + len)]
}

pred_long <- {
  n_sweeps   <- length(res$lengths)
  species_id <- as.integer(res$species_indices)
  run_index  <- ave(seq_len(n_sweeps), species_id, FUN = function(ix) seq_along(ix))
  rows <- vector("list", n_sweeps)
  for (i in seq_len(n_sweeps)) {
    n_vec <- get_sweep_n(res, i)
    if (!length(n_vec)) next
    rows[[i]] <- data.frame(
      species_id                = species_id[i],
      run                       = run_index[i],
      sample_pos_within_species = seq_along(n_vec),
      predicted_value           = as.numeric(n_vec),
      stringsAsFactors = FALSE
    )
  }
  bind_rows(rows)
}

# -------- 2) Read the existing mapping (no remapping) --------
map_df <- read_csv(species_map_csv, show_col_types = FALSE) %>%
  mutate(species_id = as.integer(species_id))

pred_named <- pred_long %>%
  left_join(map_df, by = "species_id") %>%
  rename(species = species)

# -------- 3) Load counts & ground truth; prep long tables --------
biospikeins <- c("Blattidae_cluster1","Drosophilidae_cluster1",
                 "Drosophilidae_cluster2","Drosophilidae_cluster3",
                 "Gryllidae_cluster1","Gryllidae_cluster2")

ground_truth <- read.table(gt_path, header = TRUE)
ground_truth_nobio <- ground_truth[!(ground_truth$Species %in% biospikeins), ]
gt_core <- ground_truth_nobio %>%
  transmute(species = Species, sample = Sample_no, true_value = Freq)

counts <- read.table(counts_path, sep = ";", header = TRUE, check.names = FALSE)
if (nrow(counts) >= 2) counts <- counts[-(1:2), ]
counts <- counts[!(counts$cluster %in% biospikeins), , drop = FALSE]

# Stable sort only for long-table build; mapping already fixed
df_sorted <- counts[order(counts$cluster), ]

df_long <- df_sorted %>%
  pivot_longer(cols = -cluster, names_to = "sample", values_to = "count") %>%
  mutate(sample = as.numeric(sub(".*H(\\d+).*", "\\1", sample))) %>%
  rename(species = cluster)

df_merged <- df_long %>%
  left_join(gt_core %>% select(species, sample, true_value), by = c("species","sample"))

df_merged_nonzero <- df_merged %>% filter(count != 0) %>%
  arrange(species, sample) %>%
  group_by(species) %>%
  mutate(sample_pos_within_species = dplyr::row_number()) %>%
  ungroup()

df_merge_key <- df_merged_nonzero %>%
  select(species, sample, true_value, sample_pos_within_species)

# -------- 5) Align predictions to counts per species --------
counts_k <- df_merged_nonzero %>% count(species, name = "k_nonzero_counts")

sweep_k <- pred_named %>%
  group_by(species, run) %>%
  summarise(k_pred = max(sample_pos_within_species), .groups = "drop") %>%
  left_join(counts_k, by = "species")

cat("\n--- COUNTS side ---\n")
k_counts <- df_merged_nonzero %>% count(species, name="k_nonzero_counts")
cat("Total non-zero COUNT rows (expected baseline): ",
    sum(k_counts$k_nonzero_counts), "\n")

cat("\n--- PRED side, before alignment ---\n")
pred_before <- pred_named %>% count(species, run, name="n_pred_rows")
cat("Total PRED rows before alignment: ", sum(pred_before$n_pred_rows), "\n")

cat("\n--- Alignment diffs (k_pred - k_nonzero_counts) ---\n")
align_diff <- sweep_k %>%
  mutate(diff = k_pred - k_nonzero_counts)
print(align_diff %>% count(diff, name="n_runs_by_diff") %>% arrange(diff))

cat("\nRuns dropped by STRICT (k_pred != k_nonzero_counts): ",
    align_diff %>% filter(k_pred != k_nonzero_counts) %>% nrow(), "\n")

# Alignment
if (alignment_mode == "strict") {
  good_sweeps <- sweep_k %>% filter(k_pred == k_nonzero_counts) %>% select(species, run)
  pred_aligned <- pred_named %>%
    inner_join(good_sweeps, by = c("species","run")) %>%
    left_join(counts_k, by = "species") %>%
    filter(sample_pos_within_species <= k_nonzero_counts) %>%
    select(-k_nonzero_counts)
} else if (alignment_mode == "truncate") {
  pred_aligned <- pred_named %>%
    left_join(counts_k, by = "species") %>%
    filter(sample_pos_within_species <= k_nonzero_counts) %>%
    select(-k_nonzero_counts)
} else {
  stop("alignment_mode must be 'strict' or 'truncate'")
}

# -------- 6) Sweep aggregation (collapse across runs) --------
build_sweep_meta <- function(res) {
  n_sweeps <- length(res$lengths)
  sp_id    <- as.integer(res$species_indices)
  run_idx  <- ave(seq_len(n_sweeps), sp_id, FUN = function(ix) seq_along(ix))
  data.frame(
    species_id = sp_id,
    run        = run_idx,
    normconst  = if (!is.null(res$normconst)) as.numeric(res$normconst) else NA_real_,
    stringsAsFactors = FALSE
  )
}
sweep_meta <- build_sweep_meta(res)

pred_aligned_nc <- pred_aligned %>%
  left_join(sweep_meta, by = c("species_id","run"))

# ---- Convergence diagnostics on normConst (by VARIANCE) ----
normconst_variance_threshold <- 1

nc_stats <- pred_aligned_nc %>%
  dplyr::distinct(species_id, species, run, normconst) %>%
  dplyr::group_by(species_id, species) %>%
  dplyr::summarise(
    n_runs  = dplyr::n(),
    nc_min  = suppressWarnings(min(normconst, na.rm = TRUE)),
    nc_max  = suppressWarnings(max(normconst, na.rm = TRUE)),
    nc_range= nc_max - nc_min,
    nc_sd   = suppressWarnings(sd(normconst, na.rm = TRUE)),
    nc_var  = suppressWarnings(var(normconst, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    has_any_nc = is.finite(nc_var) & !is.na(nc_var) & n_runs >= 2
  )

nonconverged <- nc_stats %>%
  dplyr::filter(has_any_nc & nc_var >= normconst_variance_threshold) %>%
  dplyr::arrange(dplyr::desc(nc_var))

cat("\n=== Convergence check (normConst VARIANCE <", normconst_variance_threshold, ") ===\n")
if (nrow(nonconverged)) {
  bad_ids <- sort(unique(nonconverged$species_id))
  cat("Species indices failing convergence (variance >= ", normconst_variance_threshold, "): ",
      paste(bad_ids, collapse = ", "), "\n", sep = "")
  print(nonconverged, n = length(nonconverged$species_id))
} else {
  cat("All species pass the normConst variance threshold.\n")
}

# Save diagnostics
diag_dir <- file.path(out_dir, "diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(nc_stats,     file.path(diag_dir, "normconst_variation_by_species.csv"))
readr::write_csv(nonconverged, file.path(diag_dir, "nonconverged_normconst_variance_ge1.csv"))

# ---- Optional: exclude non-converged species from plots ----
exclude_nonconverged <- TRUE

if (exclude_nonconverged) {
  converged_ids <- nc_stats %>%
    filter((has_any_nc & nc_var < normconst_variance_threshold) | !has_any_nc) %>%
    distinct(species_id)

  n_before <- n_distinct(pred_aligned_nc$species_id)
  pred_aligned_nc <- pred_aligned_nc %>% semi_join(converged_ids, by = "species_id")
  message("Filtered to converged species: ",
          n_distinct(pred_aligned_nc$species_id), " of ", n_before)
}

# ===================== normConst diagnostics (one species) =====================
inspect_by        <- "name"               # "name" or "id"
inspect_species   <- "Acanthosomatidae_cluster1"
inspect_species_id<- NA_integer_

.soft_w <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return(numeric(0))
  if (all(is.na(x))) return(rep(1/length(x), length(x)))
  x[is.na(x)] <- -Inf
  m <- max(x)
  w <- exp(x - m)
  if (!is.finite(sum(w)) || sum(w) == 0) rep(1/length(x), length(x)) else w / sum(w)
}

diag_tbl <- pred_aligned_nc %>%
  { if (tolower(inspect_by) == "id" && is.finite(inspect_species_id))
      dplyr::filter(., species_id == inspect_species_id)
    else dplyr::filter(., species == inspect_species)
  } %>%
  distinct(species, species_id, run, normconst) %>%
  arrange(desc(normconst)) %>%
  mutate(
    w_softmax = .soft_w(normconst),
    rank      = row_number(),
    is_MAP    = normconst == max(normconst, na.rm = TRUE)
  )

if (nrow(diag_tbl) == 0) {
  message("No rows found for the requested species. Check 'inspect_by' and species name/id.")
} else {
  cat("\n=== normConst diagnostics ===\n")
  cat("Species: ", unique(diag_tbl$species), "   (species_id: ", unique(diag_tbl$species_id), ")\n", sep = "")
  cat("Sweeps (runs) found: ", nrow(diag_tbl), "\n", sep = "")
  cat("Highest normConst (MAP) at run: ",
      paste(diag_tbl$run[diag_tbl$is_MAP], collapse = ", "),
      " (ties indicate equal MAP)\n", sep = "")
  diag_out <- diag_tbl %>%
    select(run, normconst, w_softmax, is_MAP) %>%
    arrange(desc(normconst))
  print(as.data.frame(diag_out), row.names = FALSE)
  cat("\nWeight summary (softmax over normConst):\n")
  print(summary(diag_tbl$w_softmax))
  cat("Sum of weights: ", sum(diag_tbl$w_softmax), "\n", sep = "")
}

softmax_weights <- function(x) {
  x <- as.numeric(x)
  if (!length(x)) return(numeric(0))
  if (all(is.na(x))) return(rep(1/length(x), length(x)))
  x[is.na(x)] <- -Inf
  m <- max(x)
  w <- exp(x - m)
  if (!is.finite(sum(w)) || sum(w) == 0) rep(1/length(x), length(x)) else w / sum(w)
}

stat_mode <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  tab <- sort(table(v), decreasing = TRUE)
  as.numeric(names(tab)[1])
}

message("Rows BEFORE aggregation: ", nrow(pred_aligned_nc))

if (sweep_aggregation == "all") {
  pred_collapsed <- pred_aligned_nc
} else if (sweep_aggregation == "mean") {
  pred_collapsed <- pred_aligned_nc %>%
    group_by(species_id, sample_pos_within_species) %>%
    summarise(
      species = dplyr::first(species),
      predicted_value = mean(predicted_value, na.rm = TRUE),
      .groups = "drop"
    )
    } else if (sweep_aggregation == "wmode_normconst") {
  pred_collapsed <- pred_aligned_nc %>%
    group_by(species_id, sample_pos_within_species) %>%
    mutate(
      w  = softmax_weights(normconst),
      pv = round(predicted_value)
    ) %>%
    group_by(species_id, sample_pos_within_species, pv) %>%
    summarise(
      species = dplyr::first(species),
      wsum    = sum(w, na.rm = TRUE),
      .groups = "drop_last"
    ) %>%
    slice_max(order_by = wsum, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      species_id, sample_pos_within_species,
      species,
      predicted_value = pv
    )
} else if (sweep_aggregation == "mode") {
  pred_collapsed <- pred_aligned_nc %>%
    group_by(species_id, sample_pos_within_species) %>%
    summarise(
      species = dplyr::first(species),
      predicted_value = stat_mode(round(predicted_value)),
      .groups = "drop"
    )
} else if (sweep_aggregation == "map_normconst") {
  best_runs <- pred_aligned_nc %>%
    group_by(species_id) %>%
    slice_max(order_by = normconst, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(species_id, run)
  pred_collapsed <- pred_aligned_nc %>%
    inner_join(best_runs, by = c("species_id","run")) %>%
    select(-normconst)
} else if (sweep_aggregation == "wmean_normconst") {
  pred_collapsed <- pred_aligned_nc %>%
    group_by(species_id, sample_pos_within_species) %>%
    mutate(w = softmax_weights(normconst)) %>%
    summarise(
      species = dplyr::first(species),
      predicted_value = sum(w * predicted_value),
      .groups = "drop"
    ) %>%
    mutate(predicted_value = round(predicted_value))
} else {
  stop("Unknown sweep_aggregation: ", sweep_aggregation)
}

message("Rows AFTER aggregation (", sweep_aggregation, "): ", nrow(pred_collapsed))

cat("\n--- Merge coverage diagnostics ---\n")
merge_key <- df_merge_key %>% select(species, sample_pos_within_species)
pred_keys <- pred_collapsed %>% select(species, sample_pos_within_species) %>% distinct()

missing_from_pred <- anti_join(merge_key, pred_keys,
                               by=c("species","sample_pos_within_species"))
cat("Count rows with NO matching prediction (will be dropped by inner_join): ",
    nrow(missing_from_pred), "\n")
if (nrow(missing_from_pred)) {
  print(head(missing_from_pred, 15))
}

# ========= restrict to species that actually have predictions =========
keep_species <- pred_collapsed %>% distinct(species)
df_merge_key <- df_merge_key %>%
  semi_join(keep_species, by = "species")

# -------- 7) Final merge with truth & error computation --------
df_final <- pred_collapsed %>%
  inner_join(df_merge_key, by = c("species","sample_pos_within_species")) %>%
  mutate(
    true_value      = suppressWarnings(as.numeric(true_value)),
    predicted_value = suppressWarnings(as.numeric(predicted_value)),
    error           = predicted_value - true_value
  )

message("Final merged rows: ", nrow(df_final))

# ============================================================
# ======================= PLOTTING ===========================
# ============================================================

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- 1 & 2) Per-sample violins ----------
df_filtered_le20 <- df_final %>% filter(true_value <= 20)
common_min <- min(df_filtered_le20$error, df_final$error, na.rm = TRUE)
common_max <- max(df_filtered_le20$error, df_final$error, na.rm = TRUE)
common_range <- c(common_min, common_max)

p_by_sample_all <- ggplot(df_final, aes(x = factor(sample), y = error)) +
  geom_violin(fill = "#76B7B2", alpha = 0.6) +
  theme_minimal() +
  labs(x = "Sample", y = "Error (predicted - true)",
       title = "Distribution of Errors by Samples - all") +
  coord_cartesian(ylim = common_range) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_by_sample_le20 <- ggplot(df_filtered_le20, aes(x = factor(sample), y = error)) +
  geom_violin(fill = "#59A14F", alpha = 0.6) +
  theme_minimal() +
  labs(x = "Sample", y = "Error (predicted - true)",
       title = "Distribution of Errors by Samples - true ≤ 20") +
  coord_cartesian(ylim = common_range) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "Errors_by_sample_all.jpeg"),  p_by_sample_all,  width = 12, height = 6)
ggsave(file.path(out_dir, "Errors_by_sample_le20.jpeg"), p_by_sample_le20, width = 12, height = 6)

# ---------- 3) Slab + scatter ----------
err_min <- floor(min(df_final$error, na.rm = TRUE))
err_max <- ceiling(max(df_final$error, na.rm = TRUE))

df_plot <- df_final %>%
  mutate(
    sample = suppressWarnings(as.numeric(sample)),
    error  = suppressWarnings(as.numeric(error))
  ) %>%
  mutate(
    sample_f = forcats::fct_reorder(factor(sample), error, .fun = median, .na_rm = TRUE),
    x_center = as.numeric(sample_f),
    x_left   = x_center - 0.2
  )

breaks_int_full <- seq(err_min - 0.5, err_max + 0.5, by = 1)
bin_centers <- head(breaks_int_full, -1) + 0.5

bin_index_of <- function(x, breaks) {
  idx <- findInterval(x, breaks, rightmost.closed = TRUE, all.inside = TRUE)
  pmax(1L, pmin(length(breaks) - 1L, idx))
}

df_hist <- df_plot %>%
  mutate(bin_idx = bin_index_of(error, breaks_int_full)) %>%
  group_by(x_center, bin_idx) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(x_center) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    y       = bin_centers[bin_idx],
    x_start = x_center + 0.05,
    x_end   = x_center + 0.05 + 0.65 * prop
  )

tick_df <- tidyr::expand_grid(
  x_center = sort(unique(as.numeric(df_plot$x_center))),
  y        = bin_centers
) %>%
  mutate(
    x_start = x_center + 0.55,
    x_end   = x_center + 0.70
  )

p_slab <- ggplot() +
  geom_segment(
    data = tick_df,
    aes(x = x_start, xend = x_end, y = y, yend = y),
    linewidth = 0.15, colour = "grey85", inherit.aes = FALSE
  ) +
  geom_segment(
    data = df_hist,
    aes(x = x_start, xend = x_end, y = y, yend = y),
    linewidth = 1.0, lineend = "butt",
    alpha = 0.55, colour = "#4E79A7", inherit.aes = FALSE
  ) +
  geom_point(
    data = df_plot,
    mapping = aes(x = x_left, y = error),
    size = 0.55, alpha = 0.30
  ) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey20") +
  scale_x_continuous(
    breaks = sort(unique(df_plot$x_center)),
    labels = levels(df_plot$sample_f)
  ) +
  scale_y_continuous(breaks = function(lims) sort(unique(c(pretty(lims), 0)))) +
  theme_minimal() +
  labs(x = "Sample", y = "Error (predicted - true)",
       title = "Right: empirical histogram slab (blue) · Left: overlapping scatter") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, paste0("Predict_fullset_slab_scatter_", sweep_aggregation, ".jpeg")),
       p_slab, width = 10, height = 5)

# ---------- Custom bins: singletons −5..5, plus outliers "<-5" and ">5" ----------
.custom_levels <- c("<-5", as.character(-5:5), ">5")

# Display nicer minus signs (−) for negatives; singletons stay one-line
pretty_custom_labels <- function(lvls) {
  gsub("-", "\u2212", lvls, fixed = TRUE)
}

# Map numeric error to the custom bins; NA/Inf -> NA (dropped before plotting)
custom_bin_label <- function(v) {
  if (!is.finite(v)) return(NA_character_)
  v <- as.integer(round(v))
  if (v < -5L) return("<-5")
  if (v >  5L) return(">5")
  return(as.character(v))  # singleton −5..5 (including 0)
}

custom_bin_factor <- function(x, levels = .custom_levels) {
  labs <- vapply(x, custom_bin_label, character(1))
  factor(labs, levels = levels, ordered = TRUE)
}

# ---------- 4) Overall error histograms: overlaid (two-layer) ----------
col_pred      <- "#2C7FB8"
col_pred_brdr <- "#174A7E"  # darker outline for prediction
col_prior     <- "#F16913"


alpha_pred   <- 0.90
alpha_prior  <- 0.80
width_prior  <- 0.90
width_pred   <- 0.60

set.seed(42)
N_total <- nrow(df_final)
prior_counts_overall <- prior_counts_by_species(df_final, cap = 200L, seed = 42)
prior_errors_overall <- prior_counts_overall - df_final$true_value

# Coerce to integers; keep non-finite as NA; drop before plotting
df_final$error        <- as.integer(round(df_final$error))
prior_errors_overall  <- as.integer(round(prior_errors_overall))

pred_all <- df_final %>%
  mutate(bin = custom_bin_factor(error)) %>%
  filter(!is.na(bin)) %>%
  count(bin, name = "n") %>%
  mutate(prop = n / sum(n), source = "Prediction")

prior_all <- tibble::tibble(error = prior_errors_overall) %>%
  mutate(bin = custom_bin_factor(error)) %>%
  filter(!is.na(bin)) %>%
  count(bin, name = "n") %>%
  mutate(prop = if (sum(n) > 0) n / sum(n) else 0, source = "Prior")

overall_long <- bind_rows(pred_all, prior_all)

p_all_bar_pfit <- ggplot(overall_long, aes(x = bin, y = prop, fill = source)) +
  geom_col(data = subset(overall_long, source == "Prior"),
           position = "identity", width = width_prior, alpha = alpha_prior,
           colour = NA) +
  geom_col(data = subset(overall_long, source == "Prediction"),
           position = "identity", width = width_pred,  alpha = alpha_pred,
           colour = col_pred_brdr, linewidth = 0.15) +
  scale_fill_manual(
    name = "",
    values = c("Prior" = col_prior, "Prediction" = col_pred),
    breaks = c("Prediction", "Prior")
  ) +
  guides(fill = guide_legend(override.aes = list(
    linewidth = c(0.15, 0), colour = c(col_pred_brdr, NA),
    alpha = c(alpha_pred, alpha_prior), width = c(width_pred, width_prior)
  ))) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.7),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    breaks = .custom_levels,
    labels = pretty_custom_labels,
    drop = FALSE, na.translate = FALSE
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.text.x = element_text(size = 9, lineheight = 0.95, margin = margin(t = 6))
  ) +
  labs(x = "Absolute error (predicted n - true n)",
       y = "Percent",
       title = "Error distributions across samples (overlaid histograms)")
       #subtitle = "Prior (wide, orange) under Prediction (narrow, blue)")

ggsave(file.path(out_dir, paste0("Predict_fullset_summary_", sweep_aggregation, ".jpeg")),
       p_all_bar_pfit, width = 6, height = 5)

# ---------- 5) Cutoff histogram: overlaid (two-layer) ----------
cutoff <- -20
df_cut <- df_final %>% filter(error >= cutoff)
N_cut  <- nrow(df_cut)

set.seed(43)
prior_counts_cut <- prior_counts_by_species(df_cut, cap = 200L, seed = 43)
prior_errors_cut <- prior_counts_cut - df_cut$true_value

df_cut$error       <- as.integer(round(df_cut$error))
prior_errors_cut   <- as.integer(round(prior_errors_cut))

pred_cut <- df_cut %>%
  mutate(bin = custom_bin_factor(error)) %>%
  filter(!is.na(bin)) %>%
  count(bin, name = "n") %>%
  mutate(prop = n / sum(n), source = "Prediction")

prior_cut <- tibble::tibble(error = prior_errors_cut) %>%
  mutate(bin = custom_bin_factor(error)) %>%
  filter(!is.na(bin)) %>%
  count(bin, name = "n") %>%
  mutate(prop = if (sum(n) > 0) n / sum(n) else 0, source = "Prior")

cut_long <- bind_rows(pred_cut, prior_cut)

# Counts below cutoff (use overall distributions for the text as before)
n_below_pred  <- sum(df_final$error < cutoff, na.rm = TRUE)
n_below_prior <- sum(prior_errors_overall < cutoff, na.rm = TRUE)
n_below_total <- n_below_pred + n_below_prior

# With y capped at 70%, annotate just below the cap
y_annot <- 0.67
# Cutoff is negative, so it lives in the "<-5" bin
x_cut_ix <- which(.custom_levels == "<-5")

p_all_bar_cutoff <- ggplot(cut_long, aes(x = bin, y = prop, fill = source)) +
  geom_col(data = subset(cut_long, source == "Prior"),
           position = "identity", width = width_prior, alpha = alpha_prior,
           colour = NA) +
  geom_col(data = subset(cut_long, source == "Prediction"),
           position = "identity", width = width_pred,  alpha = alpha_pred,
           colour = col_pred_brdr, linewidth = 0.15) +
  geom_vline(xintercept = x_cut_ix, linewidth = 0.6, colour = "red", linetype = "dashed") +
  annotate(
    "text",
    x = x_cut_ix + 0.25,
    y = y_annot,
    label = paste0(
      n_below_total, " outliers (< ", cutoff, ")\n",
      "Prediction: ", n_below_pred, "   ·   Prior: ", n_below_prior
    ),
    colour = "#3F3F3F",
    hjust = 0
  ) +
  scale_fill_manual(
    name = "",
    values = c("Prior" = col_prior, "Prediction" = col_pred),
    breaks = c("Prediction", "Prior")
  ) +
  guides(fill = guide_legend(override.aes = list(
    linewidth = c(0.15, 0), colour = c(col_pred_brdr, NA),
    alpha = c(alpha_pred, alpha_prior), width = c(width_pred, width_prior)
  ))) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 0.7),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    breaks = .custom_levels,
    labels = pretty_custom_labels,
    drop = FALSE, na.translate = FALSE
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.text.x = element_text(size = 9, lineheight = 0.95, margin = margin(t = 6))
  ) +
  labs(x = "Error (predicted - true) — custom bins (singletons −5..5, outliers)",
       y = "Percent",
       title = "Errors across samples (cutoff ≥ -20, overlaid histograms)",
       subtitle = "Prior (wide, orange) under Prediction (narrow, blue)")

ggsave(file.path(out_dir, paste0("Predict_fullset_summary_cutoff_", sweep_aggregation, ".jpeg")),
       p_all_bar_cutoff, width = 11, height = 4.5)

# ---------- 6 & 7) Per-species errors (violins) ----------
make_species_violin <- function(df_in, title_txt) {
  if (nrow(df_in) == 0) {
    return(ggplot() + theme_void() + labs(title = paste0(title_txt, " (no rows)")))
  }
  df_in <- df_in %>%
    mutate(species_f = forcats::fct_reorder(species, error, .fun = median, .na_rm = TRUE))
  ggplot(df_in, aes(x = species_f, y = error)) +
    geom_violin(fill = "#4E79A7", alpha = 0.55) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey20") +
    theme_minimal() +
    coord_flip() +
    labs(x = "Species", y = "Error (predicted - true)", title = title_txt) +
    theme(axis.text.y = element_text(size = 7))
}

p_species_ge5 <- make_species_violin(
  df_final %>% filter(true_value >= 5),
  paste0("Errors by species (true ≥ 5) — ", sweep_aggregation)
)

p_species_le5 <- make_species_violin(
  df_final %>% filter(true_value <= 5),
  paste0("Errors by species (true ≤ 5) — ", sweep_aggregation)
)

ggsave(file.path(out_dir, paste0("Errors_by_species_true_ge5_", sweep_aggregation, ".jpeg")),
       p_species_ge5, width = 10, height = 12)
ggsave(file.path(out_dir, paste0("Errors_by_species_true_le5_", sweep_aggregation, ".jpeg")),
       p_species_le5, width = 10, height = 12)


# -------- Terminal summary: exact error 0, 1, and -1 --------
N_total <- nrow(df_final)
n_err0  <- sum(df_final$error == 0,  na.rm = TRUE)
n_err1  <- sum(df_final$error == 1,  na.rm = TRUE)
n_errm1 <- sum(df_final$error == -1, na.rm = TRUE)

pct <- function(x) sprintf("%.1f%%", 100 * x / max(1, N_total))

cat("\n=== Exact error summary (after rounding) ===\n")
cat("Total predictions:", N_total, "\n")
cat("Error ==  0:", n_err0,  "(", pct(n_err0),  ")\n")
cat("Error ==  1:", n_err1,  "(", pct(n_err1),  ")\n")
cat("Error == -1:", n_errm1, "(", pct(n_errm1), ")\n")

df_final_pfit               <- df_final
prior_counts_overall_pfit   <- prior_counts_overall

message("Done.")
