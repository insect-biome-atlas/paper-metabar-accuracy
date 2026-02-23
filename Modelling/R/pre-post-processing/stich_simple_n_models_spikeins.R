#!/usr/bin/env Rscript
# ===============================================================
# stitch_simple_n_models.R
# For “simple” TreePPL runs that output only n (+ normConst).
#
# - Robust to different NDJSON shapes:
#     samples[[1]]$`__data__`$n
#     samples[[1]]$data$n
#     samples[[1]]$n
#     top-level: n or .n
#     nested: result$n or output$n
# - Collects per-sweep normConst when present.
# - Enforces fixed length per sweep if expected_len_per_sweep is set.
# - Explicitly opens/closes file connections.
#
# Returned object:
#   list(
#     n               = numeric (all n’s concatenated),
#     offsets         = int (start index in n for each sweep),
#     lengths         = int (length of each sweep’s n),
#     species_indices = int (species id for each sweep),
#     normconst       = numeric (normConst per sweep, NA if missing)
#   )
#
# Helpers provided:
#   sweep_count(), species_ids(), sweeps_of_species(),
#   get_sweep_n(), weights_from_normconst(), n_matrix_for_species()
# ===============================================================

suppressWarnings(suppressMessages({
  library(jsonlite)
}))

## ---------- utilities ----------

# Extract numeric values from numeric vectors or stringified vectors
.extract_nums <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.numeric(x)) return(as.numeric(x))
  s <- paste(c(x), collapse = " ")
  hits <- gregexpr("\\d+(?:\\.\\d+)?", s, perl = TRUE)
  vals <- regmatches(s, hits)[[1]]
  as.numeric(vals)
}

# Safe nested getter
.get2 <- function(x, ...) {
  cur <- x
  for (nm in list(...)) {
    if (is.null(cur)) return(NULL)
    cur <- cur[[nm]]
  }
  cur
}

# Try common locations for n
.find_n <- function(obj) {
  cands <- list(
    .get2(obj, "samples", 1, "__data__", "n"),
    .get2(obj, "samples", 1, "data", "n"),
    .get2(obj, "samples", 1, "n"),
    obj[["n"]],
    obj[[".n"]],
    .get2(obj, "result", "n"),
    .get2(obj, "output", "n")
  )
  for (c in cands) {
    v <- .extract_nums(c)
    if (length(v)) return(v)
  }
  numeric(0)
}

# Find normConst-like value (optional)
.find_normconst <- function(obj) {
  # most models use 'normConst'; some might use 'logZ'
  v <- obj$normConst
  if (is.null(v)) v <- obj$logZ
  v <- suppressWarnings(as.numeric(v))
  if (length(v) == 1 && is.finite(v)) v else NA_real_
}

# Parse one NDJSON line -> list(n=..., normConst=...) or NULL
.parse_line <- function(line_txt) {
  obj <- tryCatch(fromJSON(line_txt, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(obj)) return(NULL)

  n_vec <- .find_n(obj)
  if (!length(n_vec)) return(NULL)

  nc <- .find_normconst(obj)
  list(n = n_vec, normConst = nc)
}

## ---------- main stitcher ----------

#' Stitch NDJSON runs (one file per species), extracting n + normConst
#' @param runs_dir directory with 'pred_out_species_<id>.json' files
#' @param manifest_path JSON with field 'kept' (species ids to include)
#' @param expected_len_per_sweep integer length to keep (e.g., 15). Use NULL to allow variable.
#' @param save_rds optional path to saveRDS()
#' @param verbose logical
#' @return list(n, offsets, lengths, species_indices, normconst)
stitch_simple_n_models <- function(runs_dir,
                                   manifest_path,
                                   expected_len_per_sweep = 15,
                                   save_rds = NULL,
                                   verbose = TRUE) {
  if (!file.exists(manifest_path)) stop("manifest_path not found: ", manifest_path)
  man <- fromJSON(manifest_path)
  kept <- as.integer(man$kept)
  kept <- kept[is.finite(kept)]
  if (!length(kept)) stop("No kept species listed in manifest.")

  if (verbose) cat(sprintf("Stitching %d species from: %s\n", length(kept), runs_dir))

  files <- list.files(runs_dir, pattern = "^pred_out_species_\\d+\\.json$", full.names = TRUE)
  file_map <- setNames(files, nm = as.integer(sub("^.*pred_out_species_(\\d+)\\.json$", "\\1", files)))

  all_n <- numeric(0)
  all_offsets <- integer(0)
  all_lengths <- integer(0)
  all_species <- integer(0)
  all_nc <- numeric(0)  # normConst per sweep (NA if missing)

  for (sp in kept) {
    f <- file_map[[as.character(sp)]]
    if (is.null(f) || !file.exists(f)) {
      if (verbose) cat(sprintf("  - species %d: file not found (skipping)\n", sp))
      next
    }

    # explicit open/close for safety
    lines <- character(0)
    con <- NULL
    try({
      con <- file(f, open = "r")
      lines <- readLines(con, warn = FALSE)
    }, silent = TRUE)
    if (!is.null(con)) {
      try(close(con), silent = TRUE)
    }

    if (!length(lines)) {
      if (verbose) cat(sprintf("  - species %d: empty file (skipping)\n", sp))
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
      all_n       <- c(all_n, n_vec)
      all_offsets <- c(all_offsets, start)
      all_lengths <- c(all_lengths, length(n_vec))
      all_species <- c(all_species, sp)
      all_nc      <- c(all_nc, rec$normConst)

      sp_kept <- sp_kept + 1L
    }
    if (verbose) cat(sprintf("  - species %d: %d sweep(s) kept.\n", sp, sp_kept))
  }

  if (verbose) {
    cat("--------------------------------------------------------------\n")
    cat(sprintf("Stitched: %d sweep(s) across %d kept species.\n",
                length(all_lengths), length(unique(all_species))))
    cat(sprintf("Total n entries: %d\n", length(all_n)))
    cat(if (any(is.finite(all_nc))) "Field 'normConst' present.\n" else "Field 'normConst' absent.\n")
  }

  res <- list(
    n               = all_n,
    offsets         = as.integer(all_offsets),
    lengths         = as.integer(all_lengths),
    species_indices = as.integer(all_species),
    normconst       = as.numeric(all_nc)
  )

  if (!is.null(save_rds) && nzchar(save_rds)) {
    dir.create(dirname(save_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(res, save_rds)
    if (verbose) cat("Saved stitched RDS to: ", save_rds, "\n", sep = "")
  }
  if (verbose) cat("Done.\n")
  res
}

## ---------- helpers for downstream analysis ----------

sweep_count <- function(x) length(x$lengths)
species_ids <- function(x) sort(unique(x$species_indices))
sweeps_of_species <- function(x, sp) which(x$species_indices == sp)

get_sweep_n <- function(x, i) {
  off <- x$offsets[i]; len <- x$lengths[i]
  if (is.na(off) || is.na(len)) return(numeric(0))
  x$n[(off + 1):(off + len)]
}

weights_from_normconst <- function(x) {
  nc <- x$normconst
  if (is.null(nc) || !length(nc) || all(is.na(nc))) {
    return(rep(1/length(x$lengths), length(x$lengths)))
  }
  w <- exp(nc - max(nc, na.rm = TRUE))
  w[!is.finite(w)] <- 0
  s <- sum(w)
  if (s <= 0) rep(1/length(w), length(w)) else w/s
}

# Build a sweeps×samples matrix of n for a given species.
# If sample_count is provided, only keep sweeps of that length.
n_matrix_for_species <- function(x, sp, sample_count = NULL) {
  idx <- sweeps_of_species(x, sp)
  if (!length(idx)) return(matrix(numeric(0), nrow = 0, ncol = 0))
  if (!is.null(sample_count)) idx <- idx[x$lengths[idx] == sample_count]
  if (!length(idx)) return(matrix(numeric(0), nrow = 0, ncol = 0))

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

## ---------- quick probe (optional diagnostics) ----------

probe_lines <- function(file, n = 3) {
  cat("=== PROBE:", basename(file), "===\n")
  lines <- readLines(file, warn = FALSE, n = n)
  for (i in seq_along(lines)) {
    obj <- tryCatch(fromJSON(lines[i], simplifyVector = FALSE), error = function(e) NULL)
    cat("line", i, "\n")
    if (is.null(obj)) { cat("  <parse error>\n"); next }
    cat("  top-level keys:", paste(names(obj), collapse = ", "), "\n")
    n_vec <- .find_n(obj)
    nc    <- .find_normconst(obj)
    cat("  n length:", length(n_vec), "\n")
    cat("  normConst present:", is.finite(nc), "\n\n")
  }
  invisible(TRUE)
}

## ---------- CLI entrypoint (optional) ----------
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %in% c(3, 4)) {
    runs_dir <- args[[1]]
    manifest <- args[[2]]
    out_rds  <- args[[3]]
    expected <- if (length(args) == 4) as.integer(args[[4]]) else 15L

    stitch_simple_n_models(
      runs_dir,
      manifest,
      expected_len_per_sweep = expected,
      save_rds = out_rds
    )
  } else {
    cat("Usage: Rscript stitch_simple_n_models.R <runs_dir> <manifest.json> <out.rds> [expected_len_per_sweep]\n")
  }
}
