#!/usr/bin/env Rscript
# ================================================================
# stitch_one_species_ndjson_all_sweeps.R
# - Reads per-species NDJSON outputs from TreePPL runs
# - Each line (sweep) typically has:
#     { "samples":[{"__data__":{"m":[..], "n":[..]}, "__constructor__":"MyReturn"}],
#       "weights":[...], "normConst": <num> }
# - Stitches all species/sweeps into a compact structure:
#     list(
#       n              = numeric vector of concatenated n's,
#       offsets        = int vector (start index in n for each sweep),
#       lengths        = int vector (#elements for each sweep),
#       species_indices= int vector (species id for each sweep),
#       m              = int vector (one m per sweep, if present; else NULL),
#       normconst      = numeric vector (one per sweep; else NULL)
#     )
# - Helper functions:
#     sweep_count(x), species_ids(x), sweeps_of_species(x, sp),
#     get_sweep_n(x, i), get_sweep_m(x, i), weights_from_normconst(x),
#     n_matrix_for_species(x, sp), m_vector_for_species(x, sp)
#
# CLI usage (optional):
#   Rscript stitch_one_species_ndjson_all_sweeps.R <runs_dir> <manifest.json> <out.rds> [expected_len_per_sweep]
#
# In R:
#   source("stitch_one_species_ndjson_all_sweeps.R")
#   res <- stitch_species_autodetect(runs_dir, manifest_path, expected_len_per_sweep=15, save_rds="out.rds")
# ================================================================

suppressWarnings(suppressMessages({
  library(jsonlite)
}))

# -------- utilities ---------------------------------------------------------

# Extract numeric values even if the source is a string like "[1, 2, 3]"
.extract_nums <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.numeric(x)) return(as.numeric(x))
  s <- paste(c(x), collapse = " ")
  hits <- gregexpr("\\d+(?:\\.\\d+)?", s, perl = TRUE)
  vals <- regmatches(s, hits)[[1]]
  as.numeric(vals)
}

# Parse a single NDJSON line -> list(n= numeric(), m = integer(1)/NA, normConst = numeric(1)/NA)
.parse_line <- function(line_txt) {
  obj <- tryCatch(fromJSON(line_txt, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(obj)) return(NULL)

  n_vec <- numeric(0)
  m_val <- NA_integer_

  if (!is.null(obj$samples) && length(obj$samples) >= 1) {
    dat <- obj$samples[[1]][["__data__"]]
    if (!is.null(dat)) {
      n_vec <- .extract_nums(dat[["n"]])
      if (!is.null(dat[["m"]])) {
        m_nums <- .extract_nums(dat[["m"]])
        if (length(m_nums) >= 1) m_val <- as.integer(m_nums[1])
      }
    }
  }

  nc <- suppressWarnings(as.numeric(obj$normConst))
  if (!is.finite(nc)) nc <- NA_real_

  list(n = n_vec, m = m_val, normConst = nc)
}

# -------- stitcher ----------------------------------------------------------

#' Stitch NDJSON runs (one file per species), autodetecting m and normConst
#' @param runs_dir directory with files named 'pred_out_species_<id>.json'
#' @param manifest_path JSON manifest with field 'kept' (species ids)
#' @param expected_len_per_sweep integer or NULL; if given, keep only sweeps whose n-length == this
#' @param save_rds optional path to save the stitched list
#' @return list(n, offsets, lengths, species_indices, m, normconst)
stitch_species_autodetect <- function(runs_dir,
                                      manifest_path,
                                      expected_len_per_sweep = NULL,
                                      save_rds = NULL) {
  if (!file.exists(manifest_path)) stop("manifest_path not found: ", manifest_path)
  man <- fromJSON(manifest_path)
  kept <- as.integer(man$kept)
  kept <- kept[is.finite(kept)]
  if (!length(kept)) stop("No kept species listed in manifest.")

  cat(sprintf("Stitching %d species from: %s\n", length(kept), runs_dir))

  files <- list.files(runs_dir, pattern = "^pred_out_species_\\d+\\.json$", full.names = TRUE)
  file_map <- setNames(files, nm = as.integer(sub("^.*pred_out_species_(\\d+)\\.json$", "\\1", files)))

  all_n <- numeric(0)
  all_offsets <- integer(0)
  all_lengths <- integer(0)
  all_species <- integer(0)
  all_m <- integer(0)   # one per sweep
  all_nc <- numeric(0)  # one per sweep

  for (sp in kept) {
    f <- file_map[[as.character(sp)]]
    if (is.null(f) || !file.exists(f)) {
      cat(sprintf("  - species %d: file not found (skipping)\n", sp))
      next
    }

    lines <- readLines(f, warn = FALSE)
    if (!length(lines)) {
      cat(sprintf("  - species %d: empty file (skipping)\n", sp))
      next
    }

    sp_kept <- 0L
    for (ln in lines) {
      rec <- .parse_line(ln)
      if (is.null(rec)) next
      n_vec <- rec$n
      if (!length(n_vec)) next

      if (!is.null(expected_len_per_sweep)) {
        if (length(n_vec) != expected_len_per_sweep) next
      }

      start <- length(all_n)
      all_n <- c(all_n, n_vec)
      all_offsets <- c(all_offsets, start)
      all_lengths <- c(all_lengths, length(n_vec))
      all_species <- c(all_species, sp)

      all_m  <- c(all_m, if (is.na(rec$m)) NA_integer_ else as.integer(rec$m))
      all_nc <- c(all_nc, if (is.na(rec$normConst)) NA_real_ else rec$normConst)

      sp_kept <- sp_kept + 1L
    }
    cat(sprintf("  - species %d: %d sweep(s) kept.\n", sp, sp_kept))
  }

  cat("--------------------------------------------------------------\n")
  cat(sprintf("Stitched: %d sweep(s) across %d kept species.\n",
              length(all_lengths), length(unique(all_species))))
  cat(sprintf("Total n entries: %d\n", length(all_n)))
  cat(if (any(!is.na(all_m))) "Field 'm' present.\n" else "Field 'm' absent.\n")
  cat(if (any(!is.na(all_nc))) "Field 'normConst' present.\n" else "Field 'normConst' absent.\n")

  res <- list(
    n               = all_n,
    offsets         = as.integer(all_offsets),
    lengths         = as.integer(all_lengths),
    species_indices = as.integer(all_species),
    m               = if (any(!is.na(all_m))) as.integer(all_m) else NULL,
    normconst       = if (any(!is.na(all_nc))) as.numeric(all_nc) else NULL
  )

  if (!is.null(save_rds) && nzchar(save_rds)) {
    dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(res, save_rds)
    cat("Saved stitched RDS to: ", save_rds, "\n", sep = "")
  }
  cat("Done.\n")
  res
}

# -------- helpers (indexing & weights) --------------------------------------

# number of sweeps (lines) kept
sweep_count <- function(x) length(x$lengths)

# unique species present
species_ids <- function(x) sort(unique(x$species_indices))

# which sweep indices belong to species sp
sweeps_of_species <- function(x, sp) which(x$species_indices == sp)

# get n for sweep i as a numeric vector
get_sweep_n <- function(x, i) {
  off <- x$offsets[i]; len <- x$lengths[i]
  if (is.na(off) || is.na(len)) return(numeric(0))
  x$n[(off + 1):(off + len)]
}

# get m for sweep i (may be NA)
get_sweep_m <- function(x, i) {
  if (is.null(x$m)) return(NA_integer_)
  x$m[i]
}

# convert normconst to normalized weights per sweep
weights_from_normconst <- function(x) {
  if (is.null(x$normconst)) return(rep(1/length(x$lengths), length(x$lengths)))
  nc <- x$normconst
  if (!length(nc)) return(rep(1/length(x$lengths), length(x$lengths)))
  w <- exp(nc - max(nc, na.rm = TRUE))
  w[!is.finite(w)] <- 0
  s <- sum(w)
  if (s <= 0) rep(1/length(w), length(w)) else w/s
}

# build a sweeps×samples matrix of n for a given species
# If sample_count specified, keeps only sweeps whose length == sample_count
n_matrix_for_species <- function(x, sp, sample_count = NULL) {
  idx <- sweeps_of_species(x, sp)
  if (!length(idx)) return(matrix(numeric(0), nrow = 0, ncol = 0))
  if (!is.null(sample_count)) idx <- idx[x$lengths[idx] == sample_count]
  if (!length(idx)) return(matrix(numeric(0), nrow = 0, ncol = 0))

  # construct row-wise, enforcing consistent column count
  ncols <- x$lengths[idx[1]]
  mats <- lapply(idx, function(i) {
    v <- get_sweep_n(x, i)
    if (length(v) != ncols) return(NULL)
    v
  })
  mats <- Filter(Negate(is.null), mats)
  if (!length(mats)) return(matrix(numeric(0), nrow = 0, ncol = 0))
  do.call(rbind, mats)
}

# vector of m's for species (one per sweep)
m_vector_for_species <- function(x, sp, sample_count = NULL) {
  idx <- sweeps_of_species(x, sp)
  if (!length(idx)) return(integer(0))
  if (!is.null(sample_count)) idx <- idx[x$lengths[idx] == sample_count]
  if (!length(idx)) return(integer(0))
  if (is.null(x$m)) return(rep(NA_integer_, length(idx)))
  x$m[idx]
}

# -------- CLI entrypoint (optional) -----------------------------------------
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %in% c(3, 4)) {
    runs_dir <- args[[1]]
    manifest <- args[[2]]
    out_rds  <- args[[3]]
    expected <- if (length(args) == 4) as.integer(args[[4]]) else NULL

    stitch_species_autodetect(
      runs_dir,
      manifest,
      expected_len_per_sweep = expected,
      save_rds = out_rds
    )
  } else {
    cat("Usage: Rscript stitch_one_species_ndjson_all_sweeps.R <runs_dir> <manifest.json> <out.rds> [expected_len_per_sweep]\n")
  }
}
