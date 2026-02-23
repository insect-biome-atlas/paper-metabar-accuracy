#!/usr/bin/env Rscript
# ======================================================================
# stitch_predictions_fullset.R  — counts+map aware stitcher (index-based mapping)
#
# Assumptions:
# - species_map.csv was built index-first from the counts file with the same
#   preprocessing TreePPL used (drop first rows + spike-ins), and THEN
#   optional zero-filter AFTER mapping.
# - Therefore `species_id` == the TreePPL index (row number after preprocessing).
#
# This script:
#   1) Loads species_map.csv for names/auditing.
#   2) Loads counts and applies the SAME preprocessing to compute per-species
#      expected_k (= # non-zero positions across samples) for each species_id.
#   3) Stitches runs (pred_out_species_<id>.json[.gz]) and enforces per-species
#      expected_k using length_mode = "strict" | "truncate".
#   4) Emits diagnostics and the stitched RDS used downstream.
# ======================================================================

suppressWarnings(suppressMessages({
  library(jsonlite)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(here)
}))

repo_root  <- here::here()
pred_out   <- file.path(repo_root, "Results/step2_predictions/prediction_outputs/fullset")
stitched   <- file.path(repo_root, "Results/step2_predictions/stitched/fullset")

# --------------------------- CONFIG ----------------------------------

#inform_pfit:
runs_dir          <- file.path(pred_out, "inform_pfit/runs")
out_rds           <- file.path(stitched, "inform_pfit/stiched_pred_fullset.rds")

#inform:
#runs_dir          <- file.path(pred_out, "inform/runs")
#out_rds           <- file.path(stitched, "inform/stiched_pred_fullset.rds")

#uninform:
#runs_dir          <- file.path(pred_out, "uninform/runs")
#out_rds           <- file.path(stitched, "uninform/stiched_pred_fullset.rds")

#inform_adjust:
#runs_dir          <- file.path(pred_out, "inform_adjust/runs")
#out_rds           <- file.path(stitched, "inform_adjust/stiched_pred_fullset.rds")

#############################################

species_map_csv   <- file.path(repo_root, "Data/species_map.csv")
counts_path       <- file.path(repo_root, "Data/cleaned_nochimera_MATCHED_cluster_counts_ELA001_HOMOGEN.csv")

# enforce per-species length using counts (# non-zero positions)
length_mode              <- "strict"   # "strict" or "truncate"
auto_detect_len_global   <- TRUE       # fallback only if a run id isn't in counts
global_expected_len      <- NULL

# counts preprocessing — MUST match how the map was made / TreePPL saw the data
drop_first_n_rows <- 2
drop_spikeins     <- TRUE
biospikeins <- c("Blattidae_cluster1","Drosophilidae_cluster1",
                 "Drosophilidae_cluster2","Drosophilidae_cluster3",
                 "Gryllidae_cluster1","Gryllidae_cluster2")

# ---------------------------------------------------------------------

# ---------- utilities ----------
.extract_nums <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.numeric(x)) return(as.numeric(x))
  if (is.list(x)) {
    flat <- unlist(x, recursive = TRUE, use.names = FALSE)
    return(.extract_nums(flat))
  }
  s <- paste(c(x), collapse = " ")
  hits <- gregexpr("\\d+(?:\\.\\d+)?", s, perl = TRUE)
  vals <- regmatches(s, hits)[[1]]
  if (length(vals)) as.numeric(vals) else numeric(0)
}

.parse_record <- function(obj) {
  n_vec <- numeric(0); m_val <- NA_integer_; nc <- NA_real_
  if (!is.null(obj$samples) && length(obj$samples) >= 1) {
    dat <- obj$samples[[1]][["__data__"]]
    if (!is.null(dat)) {
      n_vec  <- .extract_nums(dat[["n"]])
      m_nums <- .extract_nums(dat[["m"]]); if (length(m_nums) >= 1) m_val <- as.integer(m_nums[1])
    }
  }
  nc_try <- suppressWarnings(as.numeric(obj$normConst))
  if (is.finite(nc_try)) nc <- nc_try
  list(n = n_vec, m = m_val, normConst = nc)
}

.read_species_file <- function(fp) {
  sz <- suppressWarnings(file.info(fp)$size)
  if (is.na(sz) || sz == 0) return(list())
  open_con <- function(path) if (grepl("\\.gz$", path, TRUE)) gzfile(path, "rt") else file(path, "rt")
  con <- open_con(fp)
  lines <- tryCatch(readLines(con, warn = FALSE), error = function(e) character(0))
  try(close(con), silent = TRUE)
  if (!length(lines)) return(list())

  # Try whole-file JSON first (object or array)
  full_txt <- paste(lines, collapse = "\n")
  try_full <- tryCatch(jsonlite::fromJSON(full_txt, simplifyVector = FALSE), error = function(e) NULL)
  if (!is.null(try_full)) {
    if (!is.list(try_full)) return(list())
    if (!is.null(try_full$normConst) || !is.null(try_full$samples)) {
      recs <- list(try_full)
    } else {
      recs <- try_full
      if (length(recs) && !is.null(names(recs)) && !is.list(recs[[1]])) recs <- as.list(recs)
    }
    return(lapply(recs, .parse_record))
  }

  # Fallback: NDJSON
  out <- vector("list", length(lines)); k <- 0L
  for (ln in lines) {
    ln_trim <- trimws(ln); if (!nzchar(ln_trim)) next
    obj <- tryCatch(jsonlite::fromJSON(ln_trim, simplifyVector = FALSE), error = function(e) NULL)
    if (is.null(obj)) next
    rec <- .parse_record(obj)
    k <- k + 1L; out[[k]] <- rec
  }
  if (k == 0L) list() else out[seq_len(k)]
}

.modal_int <- function(x) {
  x <- as.integer(x); x <- x[is.finite(x)]
  if (!length(x)) return(NA_integer_)
  tab <- table(x); as.integer(names(tab)[which.max(tab)])
}

# ---------- expected_k from counts (index-based) ----------
expected_k_from_counts <- function(counts_path, drop_first_n_rows, drop_spikeins, biospikeins) {
  counts <- read.table(counts_path, sep = ";", header = TRUE, check.names = FALSE)
  if (drop_first_n_rows > 0 && nrow(counts) >= drop_first_n_rows) {
    counts <- counts[-seq_len(drop_first_n_rows), , drop = FALSE]
  }
  stopifnot("cluster" %in% names(counts))
  if (isTRUE(drop_spikeins)) {
    counts <- counts[!(counts$cluster %in% biospikeins), , drop = FALSE]
  }
  # expected_k = # non-zero positions across samples, per row (no sorting!)
  if (ncol(counts) > 1) {
    k_nonzero <- rowSums(counts[, -1, drop = FALSE] != 0, na.rm = TRUE)
  } else k_nonzero <- integer(nrow(counts))

  tibble(
    species_id   = seq_len(nrow(counts)),     # index after preprocessing
    species      = counts$cluster,
    expected_k   = as.integer(k_nonzero)
  )
}

# ---------- stitcher ----------
stitch_species_counts_aware <- function(runs_dir,
                                        out_rds,
                                        species_map_csv,
                                        counts_path,
                                        length_mode = c("strict","truncate")[1],
                                        drop_first_n_rows = 2,
                                        drop_spikeins = TRUE,
                                        biospikeins = character(),
                                        auto_detect_len_global = TRUE,
                                        global_expected_len = NULL) {

  stopifnot(dir.exists(runs_dir))
  files <- list.files(runs_dir, pattern = "^pred_out_species_\\d+\\.json(\\.gz)?$", full.names = TRUE)
  if (!length(files)) stop("No files 'pred_out_species_<id>.json[.gz]' found in: ", runs_dir)

  run_ids <- as.integer(sub("^.*pred_out_species_(\\d+)\\.json(\\.gz)?$", "\\1", files))
  file_map <- setNames(files, nm = as.character(run_ids))

  # Load mapping (index-based) for names (may be zero-filtered)
  map_df <- readr::read_csv(species_map_csv, show_col_types = FALSE) %>%
    mutate(species_id = as.integer(species_id))

  # Compute expected_k for every index using *same preprocessing* as map/TreePPL
  exp_tbl <- expected_k_from_counts(counts_path, drop_first_n_rows, drop_spikeins, biospikeins)
  expected_k_by_id <- setNames(exp_tbl$expected_k, as.character(exp_tbl$species_id))

  # Diagnostics baseline: which IDs exist where?
  ids_expected   <- as.integer(exp_tbl$species_id)              # all indices from counts
  ids_in_map     <- sort(unique(as.integer(map_df$species_id))) # those that survived zero-filter in map
  ids_in_runs    <- sort(unique(as.integer(names(file_map))))

  ids_runs_not_in_counts <- setdiff(ids_in_runs, ids_expected)
  ids_counts_not_in_runs <- setdiff(ids_expected, ids_in_runs)
  ids_map_not_in_runs    <- setdiff(ids_in_map, ids_in_runs)

  out_dir <- dirname(out_rds)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  write_csv(tibble(species_id = ids_runs_not_in_counts),
            file.path(out_dir, "diag_runs_ids_not_in_counts.csv"))
  write_csv(tibble(species_id = ids_counts_not_in_runs) %>% left_join(exp_tbl, by = "species_id"),
            file.path(out_dir, "diag_counts_ids_missing_files.csv"))
  write_csv(tibble(species_id = ids_map_not_in_runs) %>% left_join(map_df, by = "species_id"),
            file.path(out_dir, "diag_map_ids_missing_files.csv"))

  # Global modal fallback (only used for runs IDs outside counts range)
  if (isTRUE(auto_detect_len_global) && is.null(global_expected_len)) {
    lens <- integer(0)
    for (sp_chr in names(file_map)) {
      recs <- .read_species_file(file_map[[sp_chr]])
      if (!length(recs)) next
      lens <- c(lens, vapply(recs, function(r) length(r$n), integer(1)))
    }
    gm <- .modal_int(lens[lens > 0])
    if (is.finite(gm)) {
      global_expected_len <- gm
      cat(sprintf("Global modal length (fallback): %d\n", global_expected_len))
    }
  }

  keep_or_fix <- function(n_vec, target_k, mode) {
    L <- length(n_vec)
    if (is.na(target_k) || target_k <= 0) return(list(keep = TRUE, n = n_vec, reason = "no_expected_k"))
    if (mode == "strict") {
      if (L == target_k) list(keep = TRUE, n = n_vec, reason = "match_expected_k")
      else               list(keep = FALSE, n = n_vec, reason = sprintf("drop_len_%d_ne_%d", L, target_k))
    } else { # truncate
      if (L == target_k) return(list(keep = TRUE, n = n_vec, reason = "match_expected_k"))
      if (L >  target_k) return(list(keep = TRUE, n = n_vec[seq_len(target_k)], reason = sprintf("truncate_%d_to_%d", L, target_k)))
      pad_val <- if (length(n_vec)) tail(n_vec, 1) else 0
      list(keep = TRUE, n = c(n_vec, rep(pad_val, target_k - L)), reason = sprintf("pad_%d_to_%d", L, target_k))
    }
  }

  # Stitch
  all_n <- numeric(0); all_offsets <- integer(0); all_lengths <- integer(0)
  all_species <- integer(0); all_m <- integer(0); all_nc <- numeric(0)
  diag_rows <- list()
  warned_missing_ids <- character(0)

  cat(sprintf("Stitching %d species (mode=%s)\n", length(file_map), length_mode))
  for (sp_chr in names(file_map)) {
    sp_id <- as.integer(sp_chr)
    fp    <- file_map[[sp_chr]]
    recs  <- .read_species_file(fp)

    # expected_k: from counts; fallback if run id not in counts range
    exp_k <- unname(expected_k_by_id[as.character(sp_id)])
    if (length(exp_k) == 0L || is.na(exp_k)) {
      if (!(sp_chr %in% warned_missing_ids)) {
        message(sprintf("species_id %s not present in counts-derived indices; using global fallback.", sp_chr))
        warned_missing_ids <- c(warned_missing_ids, sp_chr)
      }
      exp_k <- global_expected_len
    }

    kept_here <- 0L
    if (!length(recs)) {
      diag_rows[[length(diag_rows)+1]] <- data.frame(
        species_id = sp_id, sweep = NA_integer_, len = NA_integer_,
        decision = "no_records", reason = "empty_or_parse_fail", stringsAsFactors = FALSE
      )
      next
    }

    swe_i <- 0L
    for (rec in recs) {
      swe_i <- swe_i + 1L
      if (is.null(rec)) next
      n_vec <- rec$n; L <- length(n_vec)
      if (!L) {
        diag_rows[[length(diag_rows)+1]] <- data.frame(
          species_id = sp_id, sweep = swe_i, len = 0L,
          decision = "skip", reason = "zero_length_n", stringsAsFactors = FALSE
        )
        next
      }

      dec <- keep_or_fix(n_vec, exp_k, length_mode)
      if (!dec$keep) {
        diag_rows[[length(diag_rows)+1]] <- data.frame(
          species_id = sp_id, sweep = swe_i, len = L,
          decision = "drop", reason = dec$reason, stringsAsFactors = FALSE
        )
        next
      }

      n_use <- dec$n
      start <- length(all_n)
      all_n       <- c(all_n, n_use)
      all_offsets <- c(all_offsets, start)
      all_lengths <- c(all_lengths, length(n_use))
      all_species <- c(all_species, sp_id)
      all_m       <- c(all_m, if (is.na(rec$m)) NA_integer_ else as.integer(rec$m))
      all_nc      <- c(all_nc, if (is.na(rec$normConst)) NA_real_ else as.numeric(rec$normConst))
      kept_here   <- kept_here + 1L

      diag_rows[[length(diag_rows)+1]] <- data.frame(
        species_id = sp_id, sweep = swe_i, len = L,
        decision = "keep", reason = dec$reason, stringsAsFactors = FALSE
      )
    }
    cat(sprintf("  - species_id %s (expected_k=%s): %d sweep(s) kept.\n",
                sp_chr, ifelse(is.na(exp_k), "NA(fallback)", as.character(exp_k)), kept_here))
  }

  cat("-----------------------------------------------------------------\n")
  cat(sprintf("Stitched: %d sweep(s) across %d species.\n",
              length(all_lengths), length(unique(all_species))))
  cat(sprintf("Total n entries: %d\n", length(all_n)))
  cat(if (any(!is.na(all_m))) "Field 'm' present in some sweeps.\n" else "Field 'm' absent.\n")
  cat(if (any(!is.na(all_nc))) "Field 'normConst' present in some sweeps.\n" else "Field 'normConst' absent.\n")

  res <- list(
    n               = all_n,
    offsets         = as.integer(all_offsets),
    lengths         = as.integer(all_lengths),
    species_indices = as.integer(all_species),
    m               = if (any(!is.na(all_m))) as.integer(all_m) else NULL,
    normconst       = if (any(!is.na(all_nc))) as.numeric(all_nc) else NULL
  )

 # ---------- diagnostics ----------
diag <- dplyr::bind_rows(diag_rows)

if (nrow(diag)) {
  diag$species_id <- as.integer(diag$species_id)

  # Join counts-derived expected_k + species name, then map species name.
  diag <- diag %>%
    left_join(exp_tbl %>% select(species_id, species, expected_k), by = "species_id") %>%
    rename(species_counts = species) %>%
    left_join(map_df %>% select(species_id, species), by = "species_id", suffix = c("", "")) %>%
    rename(species_map = species) %>%
    # (Optional) a unified name if you want one column:
    mutate(species_name = dplyr::coalesce(species_map, species_counts)) %>%
    relocate(species_id, species_name, species_counts, species_map, expected_k)

  readr::write_csv(diag, file.path(out_dir, "diag_sweeps.csv"))

  species_diag <- diag %>%
    group_by(species_id, species_name, expected_k) %>%
    summarise(
      sweeps_total = n(),
      sweeps_kept  = sum(decision == "keep", na.rm = TRUE),
      sweeps_drop  = sum(decision == "drop", na.rm = TRUE),
      kept_reasons = paste(sort(unique(reason[decision == "keep"])), collapse = ";"),
      drop_reasons = paste(sort(unique(reason[decision == "drop"])), collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(desc(sweeps_drop))

  readr::write_csv(species_diag, file.path(out_dir, "diag_species_lengths.csv"))

  bad_lengths <- diag %>%
    filter(decision == "drop") %>%
    count(species_id, species_name, expected_k, len, name = "n_sweeps") %>%
    arrange(species_id, len)

  readr::write_csv(bad_lengths, file.path(out_dir, "diag_bad_lengths.csv"))
} else {
  readr::write_csv(tibble(), file.path(out_dir, "diag_sweeps.csv"))
  readr::write_csv(tibble(), file.path(out_dir, "diag_species_lengths.csv"))
  readr::write_csv(tibble(), file.path(out_dir, "diag_bad_lengths.csv"))
}


  saveRDS(res, out_rds)
  cat("Saved stitched RDS to: ", out_rds, "\n", sep = "")

  invisible(res)
}

# --------------------------- Run it ----------------------------------

res <- stitch_species_counts_aware(
  runs_dir = runs_dir,
  out_rds = out_rds,
  species_map_csv = species_map_csv,
  counts_path = counts_path,
  length_mode = length_mode,
  drop_first_n_rows = drop_first_n_rows,
  drop_spikeins = drop_spikeins,
  biospikeins = biospikeins,
  auto_detect_len_global = auto_detect_len_global,
  global_expected_len = global_expected_len
)
