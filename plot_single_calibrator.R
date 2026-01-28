# Plot the sd of log relative error for lysate and homogenate
# Swedish IBA data and 0 to 5 spike-in calibrations

# Libraries needed
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Get the coefficient of variation values
setwd("/Users/emma/git/paper-metabar-accuracy/")
source("get_0to5_calibration_error.R")


#––– Compute and plot single-spikein performance –––––––––––––––––––––––––––––––––––––––––––––––––


# Needed labels/mapping
include_cols <- c("Sh.la","Gr.bi","Gr.su","Dr.bi","Dr.se","Dr.ja")
n_vec        <- c("Sh.la"=2, "Gr.bi"=1, "Gr.su"=1, "Dr.bi"=3, "Dr.se"=1, "Dr.ja"=1)
new_labels   <- paste0(include_cols, " (n=", n_vec[include_cols], ")")

# If bio_spikes is already defined above, you can omit this re-definition:
bio_spikes   <- c(
  "Blattidae_cluster1",    # → Sh.la
  "Gryllidae_cluster1",    # → Gr.bi
  "Gryllidae_cluster2",    # → Gr.su
  "Drosophilidae_cluster1",# → Dr.bi
  "Drosophilidae_cluster2",# → Dr.se
  "Drosophilidae_cluster3" # → Dr.ja
)

### Plot Single spike-in calibration ---------------------------

# Keep only single-calibrator rows from D/E:
D1 <- D %>%
  tidyr::pivot_longer(
    cols      = dplyr::all_of(include_cols),
    names_to  = "calibrator",
    values_to = "included"
  ) %>%
  dplyr::filter(included, num_includes == 1) %>%
  dplyr::mutate(dataset = "Lysate")

E1 <- E %>%
  tidyr::pivot_longer(
    cols      = dplyr::all_of(include_cols),
    names_to  = "calibrator",
    values_to = "included"
  ) %>%
  dplyr::filter(included, num_includes == 1) %>%
  dplyr::mutate(dataset = "Homogenate")

plot_df <- dplyr::bind_rows(D1, E1) %>%
  dplyr::mutate(
    spikein    = factor(spikein,   levels = bio_spikes,   labels = new_labels),
    calibrator = factor(calibrator, levels = include_cols, labels = new_labels)
  )

# Colours: Set1 but replace yellow with black
base_set1   <- RColorBrewer::brewer.pal(n = 6, name = "Set1")
custom_cols <- base_set1; custom_cols[6] <- "#000000"
names(custom_cols) <- new_labels

### Plot efficiency by calibrator ---------------------------
# (Residual SD of target vs. the calibrator’s overall read proportion)

# Totals & proportions per spike-in, per dataset (uses LSA and HSAB already computed above)
lys_totals <- rowSums(LSA[ , -1 ])
hsa_totals <- rowSums(HSAB[ , -1 ])

prop_df_simple <- data.frame(
  spikein = rep(include_cols, 2),
  total   = c(lys_totals, hsa_totals),
  dataset = rep(c("Lysate","Homogenate"), each = length(include_cols)),
  stringsAsFactors = FALSE
) %>%
  dplyr::group_by(dataset) %>%
  dplyr::mutate(prop = total / sum(total)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    spikein = factor(spikein, levels = include_cols, labels = new_labels)
  )

# Calibrator-level proportions (match plot_df$calibrator labels)
cal_prop <- prop_df_simple %>%
  dplyr::rename(calibrator = spikein) %>%
  dplyr::select(calibrator, dataset, prop)

# Join each residual-SD point to its calibrator’s proportion
eff_all <- plot_df %>%
  dplyr::left_join(cal_prop, by = c("calibrator", "dataset"))

p_eff_by_cal <- ggplot(eff_all,
                       aes(x = prop, y = stddev,
                           colour = calibrator, shape = spikein)) +
  geom_jitter(width = 0.002, height = 0, size = 3, alpha = 0.8) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = custom_cols, drop = FALSE) +
  scale_shape_discrete(drop = FALSE) +
  labs(
   # title  = "Calibration efficiency\nresidual SD of target vs calibrator read-proportion",
    x      = "Proportion of total reads of the calibrator spike-in",
    y      = "Residual SD (log10) of the target spike-in",
    colour = "Calibrator",
    shape  = "Target"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")



ggsave("Figures/Fig_spikein_efficiency_by_calibrator.pdf",
       plot = p_eff_by_cal, width = 8, height = 4)

##################################