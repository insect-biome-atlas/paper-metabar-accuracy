# ------------------------------------------------------------
# Full script: boxplots with spike-ins left, others right
# - Calibrate by total spike-in reads per sample
# - Keep clusters in ≥10 samples (counts>0)
# - Boxes require ≥10 paired (counts∩GT) samples
# - Has homogenate and lysate sections
# ------------------------------------------------------------

# 0) packages
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# 1) Define the six biological spike-ins (cluster IDs)
bio_spikes <- c(
  "Blattidae_cluster1",
  "Gryllidae_cluster1",
  "Gryllidae_cluster2",
  "Drosophilidae_cluster1",
  "Drosophilidae_cluster2",
  "Drosophilidae_cluster3"
)

# (Optional) named mapping from cluster -> species printed on x-axis
# Include both spike-ins and any non-spike clusters you want to show
cluster <- c(
  "Blattidae_cluster1",
  "Gryllidae_cluster1",
  "Gryllidae_cluster2",
  "Drosophilidae_cluster1",
  "Drosophilidae_cluster2",
  "Drosophilidae_cluster3",
  "Entomobryidae_cluster1",
  "Muscidae_cluster28",
  "Mycetophilidae_cluster11",
  "Polleniidae_cluster1",
  "Sciaridae_cluster1",
  "Sciaridae_cluster2"
)

species <- c(
  "Shelfordella lateralis",
  "Gryllus bimaculatus",
  "Gryllodes supplicans",
  "Drosophila bicornuta",
  "Drosophila serrata",
  "Drosophila jambulina",
  "Entomobrya_X",
  "Polietes nigrolimbata",
  "Phronia basalis",
  "Pollenia vagabunda",
  "Corynoptera verrucifera",
  "Scatopsciara atomaria"
)

names_map <- data.frame(cluster, species, stringsAsFactors = FALSE)

# 2) Read & filter ground truth (used for Freq values)
ground_truth <- fread("Ground_truth_ordered.tsv", header = TRUE)
gt_filt      <- ground_truth %>% filter(Freq > 0)

# ------------------------------------------------------------
# Helper: build plot for one matrix of counts
# ------------------------------------------------------------
make_plot <- function(
  counts_path,
  sep,
  title_text,
  subtitle_text,
  out_file,
  drop_first_two_rows = FALSE
) {
  # 3) read counts
  counts_raw <- read.table(counts_path, sep = sep, header = TRUE, check.names = FALSE)
  stopifnot(colnames(counts_raw)[1] == "cluster")
  counts <- counts_raw
  if (drop_first_two_rows) {
    # (often rows 1–2 are synthetic spike-ins; keep biological ones)
    counts <- counts_raw[-c(1, 2), ]
  }

  # ------------------------------------------------------------
  # A) "occurs in >=10 samples" BY COUNTS (non-zero raw)
  # ------------------------------------------------------------
  counts_long_raw <- counts %>%
    pivot_longer(
      cols      = -cluster,
      names_to  = "Sample",
      values_to = "reads_raw"
    ) %>%
    filter(!is.na(reads_raw), reads_raw > 0)

  clusters_keep_counts <- counts_long_raw %>%
    group_by(cluster) %>%
    summarise(n_samples = n_distinct(Sample), .groups = "drop") %>%
    filter(n_samples >= 10) %>%
    pull(cluster)

  # ------------------------------------------------------------
  # B) Calibrate counts by total spike-in reads per sample
  # ------------------------------------------------------------
  spike_sums <- colSums(counts[counts$cluster %in% bio_spikes, -1, drop = FALSE], na.rm = TRUE)
  spike_sums[spike_sums == 0] <- NA_real_  # avoid dividing by 0

  counts_cal <- counts
  counts_cal[, -1] <- sweep(counts[, -1, drop = FALSE], 2, spike_sums, FUN = "/")

  # 4) long format (calibrated)
  sample_names <- names(counts_cal)[-1]
  counts_long_cal <- counts_cal %>%
    pivot_longer(
      cols      = -cluster,
      names_to  = "Sample",
      values_to = "read_counts_cal"
    ) %>%
    mutate(Sample_no = match(Sample, sample_names))

  # 5) join with ground truth and compute reads per specimen
  df_cal <- counts_long_cal %>%
    left_join(
      gt_filt %>% select(Sample_no, Species, Freq) %>% rename(cluster = Species),
      by = c("Sample_no", "cluster")
    ) %>%
    filter(!is.na(Freq), Freq > 0) %>%
    filter(!is.na(read_counts_cal)) %>%
    mutate(reads_per_specimen = read_counts_cal / Freq) %>%
    filter(is.finite(reads_per_specimen), reads_per_specimen > 0)

  # 6) keep spike-ins + non-spike clusters with >=10 samples by counts
  df_plot <- df_cal %>%
    filter(cluster %in% c(bio_spikes, clusters_keep_counts)) %>%
    mutate(
      cluster_type = ifelse(cluster %in% bio_spikes, "spike-in", "non-spike")
    )

  # Gate: >=10 paired points (same sample has counts>0 & Freq>0)
  pairs_per_cluster <- df_plot %>%
    group_by(cluster) %>%
    summarise(n_pairs = n_distinct(Sample_no), .groups = "drop") %>%
    filter(n_pairs >= 10)

  df_plot <- df_plot %>%
    inner_join(pairs_per_cluster, by = "cluster")

  # 7) ORDERING: spike-ins first (by your bio_spikes order), others alphabetically
  # Build cluster order baseline
  cluster_order <- c(
    intersect(bio_spikes, unique(df_plot$cluster)),                 # keep only present spike-ins, in your specified order
    sort(setdiff(unique(df_plot$cluster), bio_spikes))              # then non-spikes alphabetically by cluster id
  )

  # Map that cluster order to the species actually present in df_plot
  # (merge ensures we have species labels for the clusters you care about)
  df_plot <- merge(df_plot, names_map, by = "cluster")

  # Spike species (present ones, in the order of bio_spikes)
  spike_species <- names_map %>%
    filter(cluster %in% cluster_order) %>%
    mutate(order_idx = match(cluster, cluster_order)) %>%
    arrange(order_idx) %>%
    pull(species) %>%
    unique()

  # Restrict to species that are actually present in df_plot
  spike_species <- spike_species[spike_species %in% unique(df_plot$species)]

  # Set species factor levels used on x-axis: spike-ins left, others right
  # For "others", derive by cluster alphabetical order (already in cluster_order after spikes)
  species_order <- df_plot %>%
    distinct(cluster, species) %>%
    mutate(order_idx = match(cluster, cluster_order)) %>%
    arrange(order_idx) %>%
    pull(species)

  df_plot$species <- factor(df_plot$species, levels = species_order)

  # 8) plot
  p <- ggplot(df_plot, aes(x = species, y = reads_per_specimen, fill = cluster_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.3, size = 0.9) +
    scale_y_log10(
      name   = "Calibrated reads per specimen",
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = comma_format()
    ) +
    scale_fill_manual(
      name   = "Cluster type",
      values = c("spike-in" = "#E69F00", "non-spike" = "#56B4E9")  # orange vs blue
    ) +
    labs(
      x = NULL,
      title = title_text,
      subtitle = subtitle_text
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "top"
    )

  # 9) save
  ggsave(out_file, plot = p, width = 14, height = 6, dpi = 300)

  invisible(p)
}

# ------------------------------------------------------------
# Homogenate (drop first two rows if they are synthetic spike-ins)
# ------------------------------------------------------------
make_plot(
  counts_path       = "cleaned_nochimera_MATCHED_cluster_counts_ELA001_HOMOGEN.csv",
  sep               = ";",
  title_text        = "Variability of calibrated reads per specimen - homogenate",
  subtitle_text     = "Spike-ins shown first (orange); non-spikes filtered by ≥10 samples (counts-based). Boxes require ≥10 paired (counts∩GT) samples.",
  out_file          = "Figures/HOM_Reads_per_specimen_var.jpeg",
  drop_first_two_rows = TRUE
)

# ------------------------------------------------------------
# Lßysate
# ------------------------------------------------------------
make_plot(
  counts_path   = "cleaned_nochimera_MATCHED_cluster_counts_ELA001_LYS.csv",
  sep           = ",",
  title_text    = "Variability of calibrated reads per specimen - lysate",
  subtitle_text = "Spike-ins shown first (orange); non-spikes filtered by ≥10 samples (counts-based). Boxes require ≥10 paired (counts∩GT) samples.",
  out_file      = "Figures/LYS_Reads_per_specimen_var.jpeg",
  drop_first_two_rows = FALSE
)
