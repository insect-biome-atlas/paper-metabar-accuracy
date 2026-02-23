## Run the two base scripts first (they should create:
## df_final_pfit, df_final_uninform,
## prior_counts_overall_pfit, prior_counts_overall_uninform,
## p_all_bar_pfit, p_all_bar_uninform, etc.)

library(here)
repo_root <- here::here()

source(file.path(repo_root, "R/plotting/Plot_full_pred_trueprior_inform_pfit.R"))
source(file.path(repo_root, "R/plotting/Plot_full_pred_trueprior_uninform.R"))

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)
library(patchwork)  # install.packages("patchwork") if needed

#------------------------------------------------------------------
# Colors
#------------------------------------------------------------------

# Prediction (same for all panels)
col_pred      <- "#2C7FB8"
col_pred_brdr <- "#174A7E"

# Prior colors: informative vs uninformative
col_prior_inf    <- "#F16913"  # original orange
col_prior_uninf  <- "#CC79A7"  # new, clearly different (magenta-ish)

#------------------------------------------------------------------
# 1) Top row: absolute error plots (reuse p_all_bar_*)
#    + give uninformative prior a different color
#------------------------------------------------------------------

p_top_left <- p_all_bar_pfit +
  ggtitle("A) Absolute error with informative prior") +
  theme(plot.title = element_text(hjust = 0))

p_top_right <- p_all_bar_uninform +
  ggtitle("B) Absolute error with uninformative prior") +
  theme(plot.title = element_text(hjust = 0))

p_top_left <- p_top_left +
  scale_fill_manual(
    name   = "",
    values = c("Prior" = col_prior_inf, "Prediction" = col_pred),
    breaks = c("Prediction", "Prior")
  )

p_top_right <- p_top_right +
  scale_fill_manual(
    name   = "",
    values = c("Prior" = col_prior_uninf, "Prediction" = col_pred),
    breaks = c("Prediction", "Prior")
  )

#------------------------------------------------------------------
# 2) Bottom row: relative error in log factors, with:
#    - common bins & x-axis for both models
#    - percent on y
#    - prior (wide bar) behind prediction (narrow bar)
#    - **bins centered on 0.0**
#------------------------------------------------------------------

# Helper to compute log10 relative errors for one model
compute_rel_vectors <- function(df_final, prior_counts) {
  tv   <- df_final$true_value
  pred <- df_final$predicted_value

  # Prediction rel error
  mask_pred  <- is.finite(tv) & tv > 0 & is.finite(pred)
  rel_pred   <- log10(pred[mask_pred] / tv[mask_pred])

  # Prior rel error
  mask_prior <- is.finite(tv) & tv > 0 & is.finite(prior_counts)
  rel_prior  <- log10(prior_counts[mask_prior] / tv[mask_prior])

  list(pred = rel_pred, prior = rel_prior)
}

rel_inf   <- compute_rel_vectors(df_final_pfit,     prior_counts_overall_pfit)
rel_uninf <- compute_rel_vectors(df_final_uninform, prior_counts_overall_uninform)

# Global range over BOTH models & BOTH sources (pred + prior)
all_rel <- c(rel_inf$pred, rel_inf$prior, rel_uninf$pred, rel_uninf$prior)
all_rel <- all_rel[is.finite(all_rel)]

if (length(all_rel) == 0) {
  stop("No finite relative error values found.")
}

# ---- symmetric range and bins with 0 as a bin CENTER ----

max_abs <- max(abs(all_rel))
if (!is.finite(max_abs) || max_abs == 0) {
  max_abs <- 0.5  # fallback to something non-degenerate
}

# Use an ODD number of bins so one bin is centered exactly at 0
n_bins    <- 31L
bin_width <- 2 * max_abs / (n_bins - 1)

# Bin centers run from -max_abs to +max_abs, including 0
bin_centers <- seq(from = -max_abs, to = max_abs, length.out = n_bins)

# Breaks are 0.5 bin-width around centers, so 0.0 is a center, not an edge
breaks_rel <- c(bin_centers - bin_width / 2,
                tail(bin_centers, 1) + bin_width / 2)

# X axis limits for relative plots (cover all bins, no cutoffs)
rel_x_lim <- range(breaks_rel)

# Helper to build pre-binned histogram data for one (model, source)
make_rel_hist_df <- function(rel_vec, source_label, model_label) {
  tibble(rel = rel_vec) %>%
    filter(is.finite(rel)) %>%
    mutate(
      bin      = cut(rel, breaks = breaks_rel, include.lowest = TRUE, right = FALSE),
      bin_id   = as.integer(bin),
      x_center = bin_centers[bin_id]
    ) %>%
    count(bin_id, x_center, name = "n") %>%
    mutate(
      prop   = n / sum(n),
      source = source_label,
      model  = model_label
    )
}

rel_long <- bind_rows(
  make_rel_hist_df(rel_inf$pred,   "Prediction", "Informative"),
  make_rel_hist_df(rel_inf$prior,  "Prior",      "Informative"),
  make_rel_hist_df(rel_uninf$pred, "Prediction", "Uninformative"),
  make_rel_hist_df(rel_uninf$prior,"Prior",      "Uninformative")
)

rel_long_inf   <- rel_long %>% filter(model == "Informative")
rel_long_uninf <- rel_long %>% filter(model == "Uninformative")

#------------------------------------------------------------------
# Build relative plots (bottom row)
#------------------------------------------------------------------

build_rel_plot <- function(df_model, prior_color, title_txt) {
  ggplot(df_model, aes(x = x_center, y = prop, fill = source)) +
    # NEW: show 0 as dotted vertical reference line
    geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.4) +

    # Prior: wide, in the back
    geom_col(
      data    = subset(df_model, source == "Prior"),
      width   = bin_width * 0.90,
      alpha   = 0.80,
      colour  = NA,
      position = "identity"
    ) +

    # Prediction: narrower, in front, with border
    geom_col(
      data    = subset(df_model, source == "Prediction"),
      width   = bin_width * 0.60,
      alpha   = 0.90,
      colour  = col_pred_brdr,
      linewidth = 0.15,
      position = "identity"
    ) +

    scale_fill_manual(
      name   = "",
      values = c("Prior" = prior_color, "Prediction" = col_pred),
      breaks = c("Prediction", "Prior")
    ) +
    guides(fill = guide_legend(override.aes = list(
      linewidth = c(0.15, 0),
      colour    = c(col_pred_brdr, NA),
      alpha     = c(0.90, 0.80),
      width     = c(bin_width * 0.60, bin_width * 0.90)
    ))) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 0.7),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_x_continuous(
      limits = rel_x_lim,
      breaks = pretty(bin_centers, n = 6),
      labels = function(x) sprintf("%.2f", x)
    ) +
    theme_minimal() +
    theme(
      legend.position   = "top",
      legend.direction  = "horizontal",
      axis.text.x       = element_text(
        size = 9,
        angle = 45,
        hjust = 1,
        lineheight = 0.95,
        margin = margin(t = 6)
      ),
      plot.title = element_text(hjust = 0)
    ) +
    labs(
      x     = "Relative error (log10(predicted n / true n))",
      y     = "Percent",
      title = title_txt
    )
}

p_bottom_left <- build_rel_plot(
  df_model    = rel_long_inf,
  prior_color = col_prior_inf,
  title_txt   = "C) Relative error with informative prior"
)

p_bottom_right <- build_rel_plot(
  df_model    = rel_long_uninf,
  prior_color = col_prior_uninf,
  title_txt   = "D) Relative error with uninformative prior"
)

#------------------------------------------------------------------
# Output 4-panel figure
#------------------------------------------------------------------

#p_4panel <- (p_top_left + p_top_right) /
#            (p_bottom_left + p_bottom_right)

p_4panel <- (p_top_left + p_top_right) /
            plot_spacer() /
            (p_bottom_left + p_bottom_right) +
            plot_layout(heights = c(1, 0.12, 1))


ggsave(file.path(repo_root, "Results/plots/fullset/Full_pred_4panels.jpeg"),
       p_4panel, width = 11, height = 8, dpi = 300)
