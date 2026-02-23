#!/usr/bin/env Rscript

# ============================================================
# rebuild_mapping_by_index.R
# - Build mapping strictly by counts-file index (the way TreePPL indexes).
# - Make the preprocessing identical to what the model used:
#     * drop_first_n_rows   (e.g., 2)
#     * drop_spikeins       (TRUE/FALSE)
# - Then (optionally) drop zero-count rows AFTER mapping (keeps gaps in IDs).
# - Validate against runs_dir and/or stitched RDS.
# - Emit clear diagnostics for missing/excess IDs/files/rows.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(here)
})

repo_root <- here::here()

# ------------- CONFIG -------------
counts_path        <- file.path(repo_root, "Data/cleaned_nochimera_MATCHED_cluster_counts_ELA001_HOMOGEN.csv")
runs_dir           <- file.path(repo_root, "Results/step2_predictions/prediction_outputs/fullset/inform_pfit/runs")
stitched_rds_path  <- file.path(repo_root, "Results/step2_predictions/stitched/fullset/inform_pfit/stiched_pred_fullset.rds")
out_map_csv        <- file.path(repo_root, "Data/species_map.csv")
diag_dir           <- file.path(repo_root, "Data/species_map_diagnostics")

# Make these match EXACTLY what TreePPL saw when it assigned IDs:
drop_first_n_rows  <- 2          # TreePPL: did you drop 2 metadata rows? (set 0 if not)
drop_spikeins      <- TRUE       # TreePPL: did you drop spike-ins before indexing?
biospikeins <- c("Blattidae_cluster1","Drosophilidae_cluster1",
                 "Drosophilidae_cluster2","Drosophilidae_cluster3",
                 "Gryllidae_cluster1","Gryllidae_cluster2")

# After mapping, should we drop all-zero rows in the final CSV?
# (This preserves gaps in species_id if some were zero across samples.)
drop_zero_after_mapping <- TRUE

dir.create(dirname(out_map_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

# ------------- LOAD + PREPROCESS COUNTS (mirror TreePPL) -------------
counts_raw <- read.table(counts_path, sep = ";", header = TRUE, check.names = FALSE)

if (drop_first_n_rows > 0 && nrow(counts_raw) >= drop_first_n_rows) {
  counts_proc <- counts_raw[-seq_len(drop_first_n_rows), , drop = FALSE]
} else counts_proc <- counts_raw

stopifnot("cluster" %in% names(counts_proc))

if (isTRUE(drop_spikeins)) {
  counts_proc <- counts_proc[!(counts_proc$cluster %in% biospikeins), , drop = FALSE]
}

# IMPORTANT: do NOT sort here — TreePPL uses the file’s order.
species_in_order <- counts_proc$cluster
N <- length(species_in_order)

# ------------- BUILD INDEX-BASED MAP -------------
# TreePPL convention: species_id_model = row index (1..N) after the preprocessing above.
map_full <- tibble(
  species_id = seq_len(N),        # model index
  species    = species_in_order,
  canonical_index = seq_len(N)    # identical here; included for clarity
)

# ------------- OPTIONAL: DROP ZERO ROWS AFTER MAPPING -------------
if (isTRUE(drop_zero_after_mapping)) {
  if (ncol(counts_proc) > 1) {
    row_nonzero <- rowSums(counts_proc[, -1, drop = FALSE] != 0, na.rm = TRUE) > 0
  } else row_nonzero <- rep(FALSE, N)

  dropped_ids <- map_full$species_id[!row_nonzero]
  if (length(dropped_ids)) {
    message("Zero-count rows dropped AFTER mapping (IDs kept as gaps): ",
            paste(head(dropped_ids, 20), collapse = ", "),
            if (length(dropped_ids) > 20) " ..." else "")
  }
  map_final <- map_full[row_nonzero, , drop = FALSE]
} else {
  map_final <- map_full
}

# ------------- VALIDATION vs runs_dir (files) -------------
# Files are named pred_out_species_<id>.json[.gz]
if (!is.null(runs_dir) && nzchar(runs_dir) && dir.exists(runs_dir)) {
  files <- list.files(runs_dir, pattern = "^pred_out_species_\\d+\\.json(\\.gz)?$", full.names = TRUE)
  run_ids <- as.integer(sub("^.*pred_out_species_(\\d+)\\.json(\\.gz)?$", "\\1", files))
  run_ids <- sort(unique(run_ids))

  # IDs expected (from counts) vs present in runs
  expected_ids <- map_full$species_id
  ids_in_runs_not_in_map   <- setdiff(run_ids, expected_ids)
  ids_in_map_not_in_runs   <- setdiff(expected_ids, run_ids)

  write_csv(tibble(species_id = ids_in_runs_not_in_map),
            file.path(diag_dir, "runs_ids_not_in_map.csv"))
  write_csv(tibble(species_id = ids_in_map_not_in_runs) %>%
              left_join(map_full, by = "species_id"),
            file.path(diag_dir, "map_ids_missing_files.csv"))

  message(sprintf("Runs vs Map: %d ids in runs not in map; %d map ids missing files.",
                  length(ids_in_runs_not_in_map), length(ids_in_map_not_in_runs)))
} else {
  message("runs_dir not provided or missing; skipping file presence diagnostics.")
}

# ------------- VALIDATION vs stitched RDS (used IDs) -------------
if (!is.null(stitched_rds_path) && nzchar(stitched_rds_path) && file.exists(stitched_rds_path)) {
  res <- readRDS(stitched_rds_path)
  if (is.list(res) && "species_indices" %in% names(res)) {
    used_ids <- sort(unique(as.integer(res$species_indices)))
    expected_ids <- map_full$species_id

    ids_used_not_in_map <- setdiff(used_ids, expected_ids)
    ids_in_map_not_used <- setdiff(expected_ids, used_ids)

    write_csv(tibble(species_id = ids_used_not_in_map),
              file.path(diag_dir, "rds_ids_not_in_map.csv"))
    write_csv(tibble(species_id = ids_in_map_not_used) %>%
                left_join(map_full, by = "species_id"),
              file.path(diag_dir, "map_ids_unused_in_rds.csv"))

    message(sprintf("RDS vs Map: %d ids used not in map; %d map ids unused in RDS.",
                    length(ids_used_not_in_map), length(ids_in_map_not_used)))
  } else {
    warning("stitched_rds_path did not contain list$species_indices; skipping RDS diagnostics.")
  }
} else {
  message("stitched_rds_path not provided/found; skipping RDS diagnostics.")
}

# ------------- WRITE MAP -------------
# Columns: species_id (TreePPL index), species (name), canonical_index (same here)
write_csv(map_final, out_map_csv)
message("Wrote index-based mapping to: ", out_map_csv)

# Also write the full pre-zero map for auditing
write_csv(map_full, file.path(diag_dir, "species_map_pre_zero.csv"))

# ------------- QUICK SPOT CHECKS -------------
# 1) Is the map in alphabetical order? (informational)
is_alpha <- is.unsorted(map_full$species) == FALSE
message("Counts order appears alphabetical: ", is_alpha)

# 2) Show a few around 100 if present
if (N >= 110) {
  window <- map_full %>% slice(95:110)
  write_csv(window, file.path(diag_dir, "window_95_110.csv"))
  message("Wrote map window [95..110] to diagnostics for quick inspection.")
}
