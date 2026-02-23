# -------------------------------------------------------------------
# TreePPL MCMC convergence pipeline
# - Reads 3+ JSON chains
# - Extracts scalar/vector params (e.g., theta, k, c, k_bio, ...)
# - Applies burn-in + thinning
# - Computes R-hat, ESS (min across chains), Geweke (worst |z|)
# - Saves tidy CSV + plots (trace, density, R-hat summary)
# -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(coda)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- helpers --------------------------------------------------

# Turn lists-of-scalars into numeric vectors
.as_numvec <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.numeric(x)) return(as.numeric(x))

  if (is.list(x)) {
    # Try collapsing list-of-scalars to numeric
    u <- lapply(x, function(el) {
      if (is.null(el)) return(NULL)
      if (is.numeric(el)) return(as.numeric(el))
      if (is.list(el) && length(el) == 1 && is.numeric(el[[1]])) return(as.numeric(el[[1]]))
      if (is.list(el) && length(el) > 0 && is.numeric(unlist(el))) return(as.numeric(unlist(el)))
      el
    })
    if (all(vapply(u, function(v) is.numeric(v) && length(v) == 1, logical(1)))) {
      return(vapply(u, function(v) as.numeric(v)[1], numeric(1)))
    }
  }

  if (is.atomic(x)) {
    suppressWarnings(xx <- as.numeric(x))
    if (all(!is.na(xx))) return(xx)
  }
  x
}

# Pad/truncate numeric vector to length L
.pad_to <- function(v, L) {
  v <- as.numeric(v)
  if (length(v) >= L) return(v[seq_len(L)])
  c(v, rep(NA_real_, L - length(v)))
}

# ---- Extractor ----------------------------

# Normalize any recognized TreePPL shapes to a list of iterations (each a param list).
.extract_iters <- function(obj) {
  # True bulk: samples$`__data__` is an object (not a list of iterations)
  if (!is.null(obj$samples) && !is.null(obj$samples$`__data__`) && !is.list(obj$samples[[1]])) {
    return(list(obj$samples$`__data__`))  # wrap single iteration
  }

  # Common: samples is list of iterations with __data__
  if (!is.null(obj$samples) && is.list(obj$samples) && length(obj$samples) >= 1 &&
      !is.null(obj$samples[[1]]$`__data__`)) {
    return(lapply(obj$samples, function(z) z$`__data__`))
  }

  # Object is a list of length 1 holding $samples list
  if (is.list(obj) && length(obj) == 1 && !is.null(obj[[1]]$samples) && is.list(obj[[1]]$samples)) {
    s <- obj[[1]]$samples
    if (length(s) >= 1 && !is.null(s[[1]]$`__data__`)) s <- lapply(s, function(z) z$`__data__`)
    return(s)
  }

  # Top-level $chains variants
  if (!is.null(obj$chains) && is.list(obj$chains)) {
    for (ch in obj$chains) {
      if (!is.null(ch$samples) && is.list(ch$samples)) {
        if (length(ch$samples) >= 1 && !is.null(ch$samples[[1]]$`__data__`)) {
          return(lapply(ch$samples, function(z) z$`__data__`))
        }
        if (!is.null(ch$samples$`__data__`) && !is.list(ch$samples[[1]])) {
          return(list(ch$samples$`__data__`))
        }
      }
    }
  }

  stop("Could not find samples in JSON. Top-level keys: ", paste(names(obj), collapse = ", "))
}

# ---- Per-iteration flattener to wide data.frame -----------------------------

.list_of_iters_to_df <- function(s) {
  all_params <- unique(unlist(lapply(s, names)))
  if (length(all_params) == 0) stop("No parameter names found within samples.")

  out_cols <- lapply(all_params, function(p) {
    vals <- lapply(s, function(it) .as_numvec(it[[p]]))
    lens <- vapply(vals, function(v) if (is.null(v)) 0L else length(v), integer(1))
    max_len <- suppressWarnings(max(lens, na.rm = TRUE))
    if (!is.finite(max_len)) max_len <- 0L

    # NULLs -> NA/NA-vector
    vals <- lapply(vals, function(v) if (is.null(v)) (if (max_len <= 1) NA_real_ else rep(NA_real_, max_len)) else v)

    if (max_len <= 1) {
      col <- vapply(vals, function(v) as.numeric(v)[1], numeric(1))
      setNames(data.frame(col, check.names = FALSE), p)
    } else {
      mat <- t(vapply(vals, function(v) .pad_to(v, max_len), numeric(max_len)))
      dfm <- as.data.frame(mat, check.names = FALSE)
      names(dfm) <- paste0(p, "_", seq_len(max_len))
      dfm
    }
  })
  names(out_cols) <- all_params

  out <- NULL
  for (p in all_params) out <- if (is.null(out)) out_cols[[p]] else cbind(out, out_cols[[p]])
  as.data.frame(out, check.names = FALSE)
}


# ---- Chain reader -----------------------------------------------------------

read_treeppl_chain <- function(path) {
  obj <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  iters <- .extract_iters(obj)
  if (!is.list(iters) || length(iters) == 0) stop("Empty samples list.")
  .list_of_iters_to_df(iters)
}

# ---- Utilities --------------------------------------------------------------

slice_burn_thin <- function(df, burn_in, thin) {
  n <- nrow(df)
  start_idx <- if (burn_in < 1) floor(burn_in * n) + 1L else min(as.integer(burn_in) + 1L, n)
  idx <- seq.int(from = start_idx, to = n, by = thin)
  df[idx, , drop = FALSE]
}

mcmc_list_for_param <- function(chains, param) {
  mcmc.list(lapply(chains, function(df) as.mcmc(df[[param]])))
}

decide_convergence <- function(rhat, ess, geweke_z, rhat_thresh, ess_thresh) {
  pass_rhat   <- is.finite(rhat)   && (rhat <= rhat_thresh)
  pass_ess    <- is.finite(ess)    && (ess  >= ess_thresh)
  pass_geweke <- is.finite(geweke_z) && (abs(geweke_z) < 2)
  list(
    status      = if (pass_rhat && pass_ess && pass_geweke) "converged" else "needs_attention",
    pass_rhat   = pass_rhat,
    pass_ess    = pass_ess,
    pass_geweke = pass_geweke
  )
}

# ---- Main pipeline ----------------------------------------------------------

run_convergence <- function(chain_paths,
                            burn_in = 0.5,
                            thin = 1,
                            rhat_thresh = 1.05,
                            ess_thresh = 400,
                            out_dir = "mcmc_convergence") {

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  plot_dir <- file.path(out_dir, "plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  message("Reading chains...")
  raw_chains <- lapply(chain_paths, read_treeppl_chain)

  # Intersect common parameters across chains
  common_params <- Reduce(intersect, lapply(raw_chains, names))
  if (length(common_params) == 0) stop("No common parameters across chains.")
  raw_chains <- lapply(raw_chains, function(df) dplyr::select(df, all_of(common_params)))

  # Burn-in + thinning
  chains <- lapply(raw_chains, slice_burn_thin, burn_in = burn_in, thin = thin)

  # Ensure equal lengths
  lens <- vapply(chains, nrow, integer(1))
  if (length(unique(lens)) != 1) stop("Chains have unequal post burn-in/thin lengths.")

  message("Computing diagnostics...")
  results <- purrr::map_dfr(common_params, function(p) {
    ml <- mcmc_list_for_param(chains, p)

    rhat <- tryCatch({
      g <- gelman.diag(ml, autoburnin = FALSE, multivariate = FALSE)
      as.numeric(g$psrf[, "Point est."])
    }, error = function(e) NA_real_)

    ess <- tryCatch({
      es <- effectiveSize(ml)  # per-chain ESS
      as.numeric(min(es))
    }, error = function(e) NA_real_)

    gz <- tryCatch({
      zs <- vapply(seq_along(ml), function(i) {
        z <- geweke.diag(ml[[i]])
        zz <- as.numeric(z$z)
        if (length(zz) == 0) NA_real_ else zz[1]
      }, numeric(1))
      zs[which.max(abs(zs))]
    }, error = function(e) NA_real_)

    dec <- decide_convergence(rhat, ess, gz, rhat_thresh, ess_thresh)

    tibble(
      parameter  = p,
      rhat       = rhat,
      ess        = ess,
      geweke_z   = gz,
      pass_rhat  = dec$pass_rhat,
      pass_ess   = dec$pass_ess,
      pass_geweke= dec$pass_geweke,
      verdict    = dec$status
    )
  })

  # ---- Plots ----
  message("Saving plots...")
  long_df <- purrr::imap_dfr(chains, ~{
    .x %>%
      mutate(.iter = dplyr::row_number(), .chain = paste0("chain", .y)) %>%
      pivot_longer(cols = all_of(common_params), names_to = "parameter", values_to = "value")
  })

  make_trace_plot <- function(df_param, param) {
    ggplot(df_param, aes(x = .iter, y = value, color = .chain)) +
      geom_line(alpha = 0.7) +
      labs(title = paste("Trace:", param), x = "Iteration", y = param) +
      theme_bw()
  }

  make_density_plot <- function(df_param, param) {
    ggplot(df_param, aes(x = value, fill = .chain, color = .chain)) +
      geom_density(alpha = 0.2) +
      labs(title = paste("Density:", param), x = param, y = "Density") +
      theme_bw()
  }

  for (p in common_params) {
    dfp <- dplyr::filter(long_df, parameter == p)
    g1 <- make_trace_plot(dfp, p)
    g2 <- make_density_plot(dfp, p)
    ggsave(file.path(plot_dir, paste0("trace_", p, ".png")),  g1, width = 7, height = 4, dpi = 150)
    ggsave(file.path(plot_dir, paste0("density_", p, ".png")), g2, width = 6, height = 4, dpi = 150)
  }

  rhat_plot <- ggplot(results, aes(x = reorder(parameter, rhat), y = rhat, fill = verdict)) +
    geom_col() +
    geom_hline(yintercept = rhat_thresh, linetype = 2) +
    coord_flip() +
    labs(title = "R-hat summary", x = NULL, y = "R-hat") +
    theme_bw()
  ggsave(file.path(plot_dir, "rhat_summary.png"), rhat_plot, width = 7, height = 8, dpi = 150)

  # Save CSV
  write.csv(results, file.path(out_dir, "convergence_summary.csv"), row.names = FALSE)

  # Console summary
  cat("\n===== Convergence summary =====\n")
  print(results %>% arrange(verdict, rhat))
  cat("\nSaved:\n  -", file.path(out_dir, "convergence_summary.csv"),
      "\n  - plots in", plot_dir, "\n")

  invisible(results)
}
