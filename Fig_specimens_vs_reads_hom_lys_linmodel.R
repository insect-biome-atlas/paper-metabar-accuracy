# =========================
# Packages
# =========================
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(broom)

# =========================
# 1) Spike-ins and mapping
# =========================
bio_spikes <- c(
  "Blattidae_cluster1",
  "Gryllidae_cluster1",
  "Gryllidae_cluster2",
  "Drosophilidae_cluster1",
  "Drosophilidae_cluster2",
  "Drosophilidae_cluster3"
)

# ---------- cluster -> species mapping ----------
cluster_map <- c(
  "Blattidae_cluster1"       = "Shelfordella lateralis",
  "Gryllidae_cluster1"       = "Gryllus bimaculatus",
  "Gryllidae_cluster2"       = "Gryllodes supplicans",
  "Drosophilidae_cluster1"   = "Drosophila bicornuta",
  "Drosophilidae_cluster2"   = "Drosophila serrata",
  "Drosophilidae_cluster3"   = "Drosophila jambulina",
  "Entomobryidae_cluster1"   = "Entomobrya_X",
  "Muscidae_cluster28"       = "Polietes nigrolimbata",
  "Mycetophilidae_cluster11" = "Phronia basalis",
  "Polleniidae_cluster1"     = "Pollenia vagabunda",
  "Sciaridae_cluster1"       = "Corynoptera verrucifera",
  "Sciaridae_cluster2"       = "Scatopsciara atomaria"
)

# =========================
# 2) Ground truth (keep only positive Freq)
# =========================
ground_truth <- fread("Ground_truth_ordered.tsv", header = TRUE)
gt0 <- ground_truth %>% filter(Freq > 0)  # Sample_no, Species, Freq ...

# =============================================================================
# Helper: build calibrated plotting frame + per-facet stats + plot/ggsave
# =============================================================================
build_and_plot <- function(counts_file, sep, out_jpeg, y_axis_lab, plot_title, drop_first_two_rows = FALSE) {
  # ---- Read counts ----
  counts_raw <- read.table(counts_file, sep = sep, header = TRUE, check.names = FALSE)
  stopifnot(colnames(counts_raw)[1] == "cluster")
  counts <- counts_raw
  if (drop_first_two_rows) {
    # drop the first two (synthetic) rows
    counts <- counts[-c(1, 2), ]
  }

  # ---- A) PAIRED-POINT KEEPERS (>10 samples where BOTH GT>0 and counts>0) ----
  sample_names_raw <- names(counts)[-1]
  counts_long_raw <- counts %>%
    pivot_longer(
      cols      = -cluster,
      names_to  = "Sample",
      values_to = "reads_raw"
    ) %>%
    mutate(Sample_no = match(Sample, sample_names_raw)) %>%
    filter(!is.na(reads_raw), reads_raw > 0)

  # Intersection with GT (paired presence)
  paired_raw_gt <- counts_long_raw %>%
    inner_join(
      gt0 %>% select(Sample_no, Species, Freq) %>% rename(cluster = Species),
      by = c("Sample_no", "cluster")
    )

  # Keep clusters with >=10 distinct samples in this intersection
  keepers_pairs <- paired_raw_gt %>%
    group_by(cluster) %>%
    summarise(n_pairs = n_distinct(Sample_no), .groups = "drop") %>%
    filter(n_pairs >= 10) %>%
    pull(cluster)

  # Narrow GT for this dataset (keepers + spikes for calibration sanity)
  gt_filt <- gt0 %>% filter(Species %in% c(keepers_pairs, bio_spikes))

  # ---- 4) Calibration by total bio spike-ins per sample ----
  spike_sums <- colSums(counts[counts$cluster %in% bio_spikes, -1, drop = FALSE], na.rm = TRUE)
  spike_sums[spike_sums == 0] <- NA_real_  # avoid division by zero

  counts_cal <- counts
  counts_cal[, -1] <- sweep(counts[, -1, drop = FALSE], 2, spike_sums, FUN = "/")

  # ---- 5) Long (calibrated) + join to GT for plotting ----
  sample_names_cal <- names(counts_cal)[-1]
  counts_long_cal <- counts_cal %>%
    pivot_longer(
      cols      = -cluster,
      names_to  = "Sample",
      values_to = "read_counts"
    ) %>%
    mutate(Sample_no = match(Sample, sample_names_cal)) %>%
    filter(!is.na(read_counts), read_counts > 0)

  df_plot_cal <- counts_long_cal %>%
    inner_join(
      gt_filt %>% select(Sample_no, Species, Freq) %>% rename(cluster = Species),
      by = c("Sample_no", "cluster")
    ) %>%
    # keep only clusters with >10 paired points; drop spike-ins from non-spike plots
    filter(cluster %in% keepers_pairs, !cluster %in% bio_spikes)

  # ---------- species labels & ordering ----------
  df_plot_cal <- df_plot_cal %>%
    mutate(
      SpeciesLabel = ifelse(
        !is.na(cluster_map[as.character(cluster)]),
        unname(cluster_map[as.character(cluster)]),
        as.character(cluster)
      )
    )

  species_levels <- df_plot_cal %>%
    distinct(cluster, SpeciesLabel) %>%
    pull(SpeciesLabel) %>%
    unique()
  df_plot_cal$SpeciesLabel <- factor(df_plot_cal$SpeciesLabel, levels = species_levels)

  # ---- Per-facet regression stats: through origin y ~ 0 + x ----
  stats_df <- df_plot_cal %>%
    group_by(SpeciesLabel) %>%
    do({
      fit <- lm(read_counts ~ 0 + Freq, data = .)
      sm  <- summary(fit)
      ci  <- suppressMessages(confint(fit, level = 0.95))
      tibble(
        slope   = coef(fit)[[1]],
        r2      = sm$r.squared,
        p_model = pf(sm$fstatistic[1], sm$fstatistic[2], sm$fstatistic[3], lower.tail = FALSE),
        se      = sm$coefficients[1, 2],
        ci_low  = ci[1, 1],
        ci_high = ci[1, 2]
      )
    }) %>%
    ungroup() %>%
    mutate(
      label = sprintf("y = %.3g·x\nR² = %.3f\n95%% CI: [%.3g, %.3g]", slope, r2, ci_low, ci_high)
    )

  # ---- Plot with regression line & annotations ----
  p_faceted_linear <- ggplot(df_plot_cal, aes(x = Freq, y = read_counts, colour = SpeciesLabel)) +
    geom_point(alpha = 0.7, size = 2) +
    # a) show a per-facet lm through origin with SE ribbon (black)
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ 0 + x, se = TRUE, color = "black") +
    # b) also draw exact fitted line using per-facet slope (blue); keeps line visible if SE ribbon is off
    geom_abline(data = stats_df, aes(slope = slope, intercept = 0), color = "blue", linewidth = 0.6) +
    # c) per-facet labels (equation, R², CI)
    geom_text(
      data = stats_df,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.05, vjust = 1.1,
      size = 3, color = "blue"
    ) +
    scale_x_continuous("Number of specimens (per cluster)") +
    scale_y_continuous(y_axis_lab, labels = comma_format()) +
    labs(title = plot_title) +
    facet_wrap(~ SpeciesLabel, ncol = 3, scales = "free") +
    guides(colour = "none") +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  # ---- Save JPEG ----
  ggsave(out_jpeg, p_faceted_linear, width = 8, height = 6, dpi = 300)

  invisible(list(data = df_plot_cal, stats = stats_df, plot = p_faceted_linear))
}

# =============================================================================
# HOMOGENATE
# =============================================================================
hom_res <- build_and_plot(
  counts_file       = "cleaned_nochimera_MATCHED_cluster_counts_ELA001_HOMOGEN.csv",
  sep               = ";",
  out_jpeg          = "Figures/HOM_Per_cluster_faceted_linear_calibrated.jpeg",
  y_axis_lab        = "Calibrated read counts - homogenate",
  plot_title        = "Homogenate Per-cluster read counts vs specimens (calibrated, linear scales)",
  drop_first_two_rows = TRUE   # drop synthetic rows 1–2 as in your original
)

# =============================================================================
# LYSATE
# =============================================================================
lys_res <- build_and_plot(
  counts_file       = "cleaned_nochimera_MATCHED_cluster_counts_ELA001_Lys.csv",
  sep               = ",",
  out_jpeg          = "Figures/LYS_Per_cluster_faceted_linear_calibrated.jpeg",
  y_axis_lab        = "Calibrated read counts - lysate",
  plot_title        = "Lysate Per-cluster read counts vs specimens (calibrated, linear scales)",
  drop_first_two_rows = FALSE  # no dropping for LYS in your original
)
