# 0) load all needed packages
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# 1) define the six bio‐spike clusters
bio_spikes <- c(
  "Blattidae_cluster1",
  "Gryllidae_cluster1",
  "Gryllidae_cluster2",
  "Drosophilidae_cluster1",
  "Drosophilidae_cluster2",
  "Drosophilidae_cluster3"
)

## 2) read & filter ground_truth (keep only positive Freq)
ground_truth <- fread("Ground_truth_ordered.tsv", header = TRUE)
gt_filt      <- ground_truth %>% filter(Freq > 0)

# 2a) keep only clusters occurring in ≥10 distinct samples
clusters_keep <- gt_filt %>%
  group_by(Species) %>%                         # assumes column is named "Species"
  summarise(n_samples = n_distinct(Sample_no), .groups = "drop") %>%
  filter(n_samples >= 10) %>%
  pull(Species)

# 2b) drop all other clusters from gt_filt
gt_filt <- gt_filt %>% filter(Species %in% clusters_keep)

# 3) total specimens per sample (for kept clusters; kept for sanity checks if needed)
dataset <- gt_filt %>%
  group_by(Sample_no) %>%
  summarise(total_Freq = sum(Freq), .groups = "drop")

# 4) read counts and drop rows 1–2 (per your note)
counts_raw <- read.table(
  "cleaned_nochimera_MATCHED_cluster_counts_ELA001_HOMOGEN.csv",
  sep = ";", header = TRUE
)
counts <- counts_raw[-c(1, 2), ]

# 4a) filter counts to the same kept clusters
counts <- counts %>% filter(cluster %in% clusters_keep)

# 5) compute calibration factors = sum of reads of the 6 bio‐spikes per sample
spike_sums <- colSums(counts[counts$cluster %in% bio_spikes, -1, drop = FALSE], na.rm = TRUE)

# 6) calibrate: divide every cluster count by its sample’s spike‐sum
counts_cal <- counts
counts_cal[, -1] <- sweep(counts[, -1, drop = FALSE], 2, spike_sums, FUN = "/")

# 10) long format + sample indices + rename to "cluster"
sample_names <- names(counts_cal)[-1]
counts_long_cal <- counts_cal %>%
  pivot_longer(
    cols      = -cluster,
    names_to  = "Sample",
    values_to = "read_counts"
  ) %>%
  mutate(
    Sample_no = match(Sample, sample_names),  # assumes sample order matches 1..N
    cluster   = cluster
  ) %>%
  filter(read_counts > 0)

# --- Build plotting DF by joining calibrated reads with ground-truth Freq ---
# gt_filt has columns: Sample_no, Species, Freq  (Species==cluster here)
df_plot_cal <- counts_long_cal %>%
  left_join(
    gt_filt %>% select(Sample_no, Species, Freq) %>% rename(cluster = Species),
    by = c("Sample_no", "cluster")
  ) %>%
  filter(!is.na(Freq))   # keep pairs actually seen in ground truth

# A) drop spike‐ins first
df_plot_cal2 <- df_plot_cal %>%
  filter(!cluster %in% bio_spikes)

# B) keep only clusters with ≥10 positive points (after dropping spike-ins)
clusters_keep2 <- df_plot_cal2 %>%
  group_by(cluster) %>%
  summarise(n_pos = n_distinct(Sample_no), .groups = "drop") %>%
  filter(n_pos >= 10) %>%
  pull(cluster)

# C) subset for plotting
df_plot_cal2 <- df_plot_cal2 %>% filter(cluster %in% clusters_keep2)

# ——————————————————————————————————————————————
# 1) Linear‐scale facetted plot
#    Each facet’s dots are a distinct color by mapping colour=cluster (one cluster per facet)
p_faceted_linear <- ggplot(df_plot_cal2, aes(x = Freq, y = read_counts, colour = cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_x_continuous("Number of specimens (per cluster)") +
  scale_y_continuous("Calibrated read counts", labels = comma_format()) +
  labs(title = "Per-cluster read counts vs specimens (calibrated, linear scales)") +
  facet_wrap(~ cluster, ncol = 3, scales = "free") +
  guides(colour = "none") +      # hide legend; color is redundant with facet
  theme_minimal(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("Per_cluster_faceted_linear.pdf", plot = p_faceted_linear, width = 9, height = 6)

# ——————————————————————————————————————————————
# 2) Log–log‐scale facetted plot
p_faceted_loglog <- ggplot(df_plot_cal2, aes(x = Freq, y = read_counts, colour = cluster)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_x_log10(
    name   = "Frequency (per cluster)",
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_log10(
    name   = "Calibrated read counts",
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = comma_format()
  ) +
  annotation_logticks(sides = "lb") +
  labs(title = "Per-cluster read counts vs frequency (calibrated, log–log scales)") +
  facet_wrap(~ cluster, ncol = 3, scales = "free") +
  guides(colour = "none") +
  theme_minimal(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("Per_cluster_faceted_loglog.pdf", plot = p_faceted_loglog, width = 9, height = 6)