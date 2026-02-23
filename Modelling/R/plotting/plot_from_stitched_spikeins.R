#!/usr/bin/env Rscript
# ===============================================================
# plot_from_stitched_multi.R  (MERGED VERSION - Best of both)
#
# Combines:
# - All core functions from original
# - Improved final plot logic from modified version
#
# Aggregation modes for combining runs (sweeps):
#   --agg all              # equal-weight average across runs
#   --agg mean             # alias of 'all'
#   --agg map              # single best run by normconst
#   --agg wmean_normconst  # softmax(normconst) weights (default)
# ===============================================================

suppressWarnings(suppressMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(gridExtra)
}))

# ---------- utils ----------
softmax <- function(x) {
  if (length(x) == 0) return(numeric(0))
  if (all(is.na(x))) return(rep(1/length(x), length(x)))
  z <- x - max(x, na.rm = TRUE)
  p <- exp(z); p / sum(p)
}

read_names <- function(src) {
  if (is.null(src) || !nzchar(src)) return(NULL)
  if (file.exists(src)) {
    x <- readLines(src, warn = FALSE)
    x <- trimws(x)
    x[nzchar(x)]
  } else {
    x <- unlist(strsplit(src, ","))
    x <- trimws(x)
    x[nzchar(x)]
  }
}

read_true_n <- function(src) {
  if (is.null(src) || !nzchar(src)) return(NULL)
  if (file.exists(src)) {
    v <- scan(src, what = numeric(), quiet = TRUE, sep = "\n")
    v[is.finite(v)]
  } else {
    nums <- unlist(strsplit(src, "[,;\\s]+"))
    as.numeric(nums[nzchar(nums)])
  }
}

make_prior_df <- function(x_vals = 1:12, n_samp = 10000L) {
  set.seed(1)
  u <- runif(n_samp)
  # two-component prior similar to earlier plots
  p <- ifelse(u < 0.65, 1, rnorm(n_samp, 0.5, 0.1))
  p[p < 0.001] <- 0.001; p[p > 0.999] <- 0.999
  data.frame(
    value = x_vals,
    prob  = sapply(x_vals, function(k) mean(dgeom(k - 1, prob = p)))
  )
}

# Extract per-species long DataFrame (all sweeps for that species)
collect_species_long <- function(res, species_id) {
  idx <- data.frame(
    sweep   = seq_along(res$lengths),
    species = res$species_indices,
    start   = res$offsets + 1L,
    end     = res$offsets + res$lengths,
    len     = res$lengths,
    m       = if (!is.null(res$m)) res$m else NA_integer_,
    wlog    = if (!is.null(res$normconst)) res$normconst else NA_real_
  )
  sub <- idx[idx$species == species_id, , drop = FALSE]
  if (!nrow(sub)) return(list(df_long = tibble(), sub = sub))

  n_list <- lapply(seq_len(nrow(sub)), \(i) res$n[sub$start[i]:sub$end[i]])
  maxL <- max(sub$len)
  mat <- t(sapply(n_list, function(v) c(v, rep(NA_real_, maxL - length(v)))))
  colnames(mat) <- paste0("V", seq_len(maxL))
  df <- as.data.frame(mat)
  df$run <- seq_len(nrow(df))

  df_long <- df |>
    tidyr::pivot_longer(dplyr::starts_with("V"), names_to = "experiment", values_to = "value") |>
    dplyr::filter(!is.na(value)) |>
    dplyr::mutate(SampleID = as.integer(stringr::str_remove(experiment, "V"))) |>
    dplyr::select(run, SampleID, value)

  list(df_long = df_long, sub = sub)
}

# Choose weights per aggregation mode
compute_weights <- function(wlog, mode = "wmean_normconst") {
  R <- length(wlog)
  mode <- tolower(mode)
  if (R == 0) return(numeric(0))
  if (mode %in% c("all","mean")) {
    rep(1/R, R)
  } else if (mode == "map") {
    if (all(is.na(wlog))) {
      w <- rep(0, R); w[1] <- 1; w
    } else {
      j <- which.max(wlog)
      w <- rep(0, R); w[j] <- 1; w
    }
  } else if (mode == "wmean_normconst") {
    if (all(is.na(wlog))) rep(1/R, R) else softmax(wlog)
  } else {
    stop("Unknown --agg mode: ", mode)
  }
}

agg_label_pretty <- function(mode) {
  c(
    all = "Equal-weight average (all runs)",
    mean = "Equal-weight average (all runs)",
    map = "MAP (best normConst run only)",
    wmean_normconst = "Weighted by normConst (softmax)"
  )[tolower(mode)]
}

# Build posterior over n with selected aggregation mode
aggregate_posterior <- function(df_long, wlog, mode, x_max) {
  if (!nrow(df_long)) return(tibble())
  w <- compute_weights(wlog, mode)
  dfw <- df_long |> left_join(
    tibble(run = seq_along(w), weight = w), by = "run"
  )
  df_post <- dfw |>
    group_by(SampleID, value) |>
    summarise(prob = sum(weight), .groups = "drop") |>
    group_by(SampleID) |>
    mutate(prob = prob / sum(prob)) |>
    ungroup()
  # fold any values > x_max into x_max for plotting axis stability
  df_post$value <- pmin(df_post$value, x_max)
  df_post |>
    group_by(SampleID, value) |>
    summarise(prob = sum(prob), .groups = "drop")
}

# Optionally, get chosen m draws for match-* models using the same weights
aggregate_m_draws <- function(sub_idx, wlog, mode, model_label) {
  if (!nrow(sub_idx) || all(is.na(sub_idx$m))) return(tibble())
  w <- compute_weights(wlog, mode)
  tibble(run = seq_len(nrow(sub_idx)),
         value = sub_idx$m,
         weights = w) |>
    filter(!is.na(value)) |>
    mutate(model = model_label)
}

# Plot n posteriors (facetted by sample)
plot_n_facets <- function(posterior_df, prior_df, true_n, spikein_name, filename, title, x_max = 12) {
  p <- ggplot(posterior_df, aes(x = value, y = prob)) +
    geom_col(fill = "grey70", width = 1, alpha = 0.9, color = "black") +
    geom_col(data = prior_df, aes(x = value, y = prob),
             width = 1, fill = NA, color = "blue", linetype = "dashed", inherit.aes = FALSE) +
    geom_vline(xintercept = true_n, color = "red", linewidth = 1) +
    facet_wrap(~ SampleID, scales = "free_y") +
    scale_x_continuous(breaks = 1:x_max, limits = c(NA, x_max)) +
    labs(title = title,
         x = "Predicted n", y = "Probability",
         caption = sprintf("Blue dashed: Prior | Red: True (n=%s) — %s", true_n, spikein_name)) +
    theme_minimal(base_size = 11) + theme(legend.position = "none")
  ggsave(filename, plot = p, width = 3000/300, height = 1550/300, dpi = 300)
}

# One model → all species plots + summaries (uses aggregation toggle)
plot_model_from_res <- function(res, model_label, outdir_model,
                                spikein_names, true_n_vec, x_max = 12,
                                is_match = FALSE,
                                agg_mode = "wmean_normconst",
                                write_weights_debug = TRUE) {

  dir.create(outdir_model, recursive = TRUE, showWarnings = FALSE)
  prior_df <- make_prior_df(1:x_max)

  samples_long_all <- dplyr::tibble()
  m_draws_all      <- dplyr::tibble()

  nspecies <- length(spikein_names)
  for (i in seq_len(nspecies)) {
    spikein_name <- spikein_names[i]
    true_n <- if (i <= length(true_n_vec) && is.finite(true_n_vec[i])) true_n_vec[i] else NA_real_

    col <- collect_species_long(res, i)
    if (!nrow(col$df_long)) next

    # Aggregate with chosen mode
    posterior_df <- aggregate_posterior(col$df_long, col$sub$wlog, agg_mode, x_max) |>
      dplyr::mutate(error = value - true_n,
                    spikein = spikein_name,
                    model = model_label)

    # per-species n plot
    fn <- file.path(outdir_model, sprintf("pred_%s_%s_spikein_%d.jpeg",
                                          gsub("\\s+","_",tolower(model_label)),
                                          tolower(agg_mode), i))
    plot_n_facets(
      posterior_df, prior_df, true_n, spikein_name, fn,
      sprintf("Predictions of n — %s\nAggregation: %s",
              model_label, agg_label_pretty(agg_mode)),
      x_max
    )

    samples_long_all <- dplyr::bind_rows(samples_long_all, posterior_df)

    # m histogram for match models (weighted by the same scheme)
    if (is_match) {
      md <- aggregate_m_draws(col$sub, col$sub$wlog, agg_mode, model_label)
      if (nrow(md)) {
        p_m <- ggplot(md, aes(x = factor(value, levels = seq_along(spikein_names), labels = spikein_names),
                              weight = weights)) +
          geom_bar(color = "black", fill = "skyblue") +
          scale_x_discrete(drop = FALSE) +
          labs(x = "Chosen spike-in (m)", y = "Weighted count",
               title = sprintf("Chosen m for %s (%s)\nAggregation: %s",
                               spikein_name, model_label, agg_label_pretty(agg_mode))) +
          theme_minimal() + theme(axis.text.x = element_text(angle = 55, hjust = 1))
        fnm <- file.path(outdir_model, sprintf("pred_%s_m_%s_spikein_%d.jpeg",
                                               gsub("\\s+","_",tolower(model_label)),
                                               tolower(agg_mode), i))
        ggsave(fnm, plot = p_m, width = 3000/300, height = 1550/300, dpi = 300)
      }
      if (nrow(md)) {
        md$spikein <- spikein_name
        m_draws_all <- dplyr::bind_rows(m_draws_all, md)
      }
    }

    # Optional: write a quick weights debug CSV once per model (for species 1)
    if (write_weights_debug && i == 1 && nrow(col$sub)) {
      wdbg <- tibble(run = seq_len(nrow(col$sub)),
                     normconst = col$sub$wlog,
                     w_all = compute_weights(col$sub$wlog, "all"),
                     w_map = compute_weights(col$sub$wlog, "map"),
                     w_wmean = compute_weights(col$sub$wlog, "wmean_normconst"))
      write.csv(wdbg, file.path(outdir_model, sprintf("weights_debug_%s.csv", tolower(agg_mode))), row.names = FALSE)
      message("Wrote weights debug to: ", file.path(outdir_model, sprintf("weights_debug_%s.csv", tolower(agg_mode))))
    }
  }

  # per-sample summary violins (errors under chosen aggregation) — WEIGHTED by prob + COLORED by spikein
  if (nrow(samples_long_all)) {
    # Facetted by sample
    p1 <- ggplot(samples_long_all, aes(x = spikein, y = error, fill = spikein)) +
      geom_violin(aes(weight = prob), trim = TRUE, alpha = 0.6, color = "black") +
      geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.8) +
      facet_wrap(~ SampleID, scales = "free_y") +
      labs(title = sprintf("Error Distribution per Sample — %s\nAggregation: %s",
                           model_label, agg_label_pretty(agg_mode)),
           x = "Spikein", y = "Error (predicted - true)") +
      theme_minimal() +
      scale_y_continuous(limits = c(-4, 16), breaks = seq(-3, 16, by = 4)) +
      theme(axis.text.x = element_text(angle = 55, hjust = 1),
            legend.position = "none")
    ggsave(file.path(outdir_model, sprintf("Error_Violin_AllSpikeins_perSample_%s_%s.jpeg",
                                           gsub("\\s+","_",tolower(model_label)),
                                           tolower(agg_mode))),
           plot = p1, width = 12, height = 8, dpi = 300)

    # Overall violin
    p2 <- ggplot(samples_long_all, aes(x = spikein, y = error, fill = spikein)) +
      geom_violin(aes(weight = prob), trim = TRUE, alpha = 0.6, color = "black") +
      geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.8) +
      labs(title = sprintf("%s — Aggregation: %s",
                           model_label, agg_label_pretty(agg_mode)),
           x = "Biological Spikein", y = "Error (predicted - true)") +
      scale_y_continuous(limits = c(-3, 12), breaks = seq(-2, 11, by = 1), minor_breaks = NULL) +
      geom_hline(yintercept = 0, linewidth = 0.25, colour = "dark grey") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "none")
    ggsave(file.path(outdir_model, sprintf("Error_Violin_By_Spikein_%s_%s.jpeg",
                                           gsub("\\s+","_",tolower(model_label)),
                                           tolower(agg_mode))),
           plot = p2, width = 8, height = 6, dpi = 300)
  }

  # combined m histogram across spikeins (match only)
  if (is_match && nrow(m_draws_all)) {
    pm <- ggplot(m_draws_all, aes(x = value, weight = weights)) +
      geom_histogram(breaks = seq(0.5, length(spikein_names) + 0.5, by = 1),
                     closed = "right", colour = "black", fill = "steelblue") +
      scale_x_continuous(breaks = 1:length(spikein_names), labels = spikein_names) +
      labs(title = sprintf("Chosen m — %s\nAggregation: %s",
                           model_label, agg_label_pretty(agg_mode)),
           x = "Spike-in used", y = "Weighted count") +
      theme_minimal() + theme(axis.text.x = element_text(angle = 55, hjust = 1))
    ggsave(file.path(outdir_model, sprintf("Hist_Chosen_m_%s_%s.jpeg",
                                           gsub("\\s+","_",tolower(model_label)),
                                           tolower(agg_mode))),
           plot = pm, width = 10, height = 4, dpi = 300)
  }

  invisible(list(samples_long = samples_long_all, m_draws = m_draws_all))
}

# ---------- IMPROVED main driver (from modified version) ----------
run_plotter_multi <- function(models, output_dir,
                              spikein_names_src = NULL, true_n_src = NULL,
                              default_names = c("Shelfordella lateralis (2)",
                                                "Drosophila bicornuta (3)",
                                                "Drosophila serrata (1)",
                                                "Drosophila jambulina (1)",
                                                "Gryllus bimaculatus (1)",
                                                "Gryllodes sigillatus (1)"),
                              default_true_n = c(2,3,1,1,1,1),
                              x_max = 12,
                              agg_mode = "wmean_normconst") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  spikein_names <- read_names(spikein_names_src)
  if (is.null(spikein_names)) spikein_names <- default_names

  true_n_vec <- read_true_n(true_n_src)
  if (is.null(true_n_vec)) true_n_vec <- default_true_n

  # align lengths
  if (length(true_n_vec) < length(spikein_names)) {
    true_n_vec <- c(true_n_vec, rep(NA_real_, length(spikein_names) - length(true_n_vec)))
  }
  if (length(spikein_names) < length(true_n_vec)) {
    spikein_names <- c(spikein_names, paste0("Spikein_", seq_len(length(true_n_vec) - length(spikein_names))))
  }

  model_order <- c("bio","art","combined","match6k","match8k")
  pretty_name <- c(bio="Biological model",
                   art="Synthetic model",
                   combined="Combined model",
                   match6k="Match model 6k",
                   match8k="Match model 8k")

  results_list <- list()

  for (key in model_order) {
    rds <- models[[key]]
    if (is.null(rds) || !nzchar(rds)) next
    if (!file.exists(rds)) {
      message(sprintf("WARN: '%s' path not found: %s (skipping)", key, rds))
      next
    }
    res <- readRDS(rds)
    is_match <- !is.null(res$m)  # match outputs include m

    outdir_model <- file.path(output_dir, key)
    results_list[[key]] <- plot_model_from_res(
      res, model_label = pretty_name[[key]], outdir_model = outdir_model,
      spikein_names = spikein_names, true_n_vec = true_n_vec,
      x_max = x_max, is_match = is_match, agg_mode = agg_mode
    )
  }

## ---------- IMPROVED FINAL VIOLIN COMBINATION (2x2 grid, drop Synthetic model) ----------
res_nonnull <- Filter(function(x) {
  !is.null(x) && !is.null(x$samples_long) && nrow(x$samples_long) > 0
}, results_list)

if (length(res_nonnull) > 0) {
  all_samples <- dplyr::bind_rows(lapply(res_nonnull, `[[`, "samples_long")) |>
    dplyr::filter(model != "Synthetic model") |>
    dplyr::filter(is.finite(error), !is.na(prob), prob >= 0)

  if (nrow(all_samples) > 0) {
    # Use only models that actually have data (prevents empty facets that can crash some ggplot2 builds)
    model_priority <- c("Biological model","Combined model","Match model 6k","Match model 8k")
    models_present <- intersect(model_priority, unique(as.character(all_samples$model)))

    # Prepend A), B), C), D) labels to subplot titles
    subplot_labels <- setNames(
      paste0(LETTERS[seq_along(models_present)], ") ", models_present),
      models_present
    )
    all_samples$model <- factor(
      subplot_labels[as.character(all_samples$model)],
      levels = unname(subplot_labels)
    )

    # Freeze spike-in order for consistent x arrangement
    spike_levels <- unique(all_samples$spikein)
    all_samples$spikein <- factor(all_samples$spikein, levels = spike_levels)

    # Boxplots only where group size >= 2
    box_df <- all_samples |>
      dplyr::group_by(model, spikein) |>
      dplyr::filter(dplyr::n() >= 2) |>
      dplyr::ungroup()
    if (nrow(box_df) > 0) {
      box_df$spikein <- factor(box_df$spikein, levels = spike_levels)
      box_layer <- ggplot2::geom_boxplot(
        data = box_df, width = 0.12, outlier.shape = NA,
        fill = "white", alpha = 0.9
      )
    } else {
      box_layer <- NULL
    }

    p_final <- ggplot(all_samples, aes(x = spikein, y = error, fill = spikein)) +
      geom_violin(aes(weight = prob), trim = TRUE, alpha = 0.65, color = "black") +
      box_layer +
      geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey40") +
      scale_y_continuous(limits = c(-3, 12), breaks = seq(-2, 11, by = 1), minor_breaks = NULL) +
      facet_wrap(~ model, ncol = 2, scales = "free_y") +   # 2x2 grid using only present models
      labs(
        title = paste0("Error (predicted − true) by Spike-in, per Model — Aggregation: ",
                       agg_label_pretty(agg_mode)),
        x = "Spike-in",
        y = "Error (predicted − true)"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 55, hjust = 1),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 14),
        plot.title  = element_text(face = "bold", size = 16)
      )

    ggsave(
      file.path(output_dir, paste0("FINAL_VIOLIN_BySpikein_ByModel_2x2_", tolower(agg_mode), ".jpeg")),
      plot = p_final, width = 12, height = 10, dpi = 300
    )
  }
}

  invisible(results_list)
}