## ------------------------------------------------------------
## Priors vs observed Geometric counts (histograms + p-panels)
## One p per Species, one k per observation
## ------------------------------------------------------------

set.seed(123)

library(here)
repo_root <- here::here()

library(MASS)      # for fitdistr (load before dplyr to avoid masking)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## 1. Read and clean observed data -----------------------------

ground_truth <- read.table(file.path(repo_root, "Data/Ground_truth_adjustedorder.tsv"), header = TRUE)

# Keep only valid geometric support
gt <- ground_truth %>%
  filter(Freq > 0)

# Confirm Species is present
stopifnot("Species" %in% names(gt))

cat("Counts of observed k:\n")
print(gt %>% count(Freq) %>% arrange(Freq))

## 1b. Read truth_vector and estimate p from data --------------

truth_vector <- read.table(file.path(repo_root, "Data/truth_vector.csv"), sep = ",")

# Coerce your input once, safely
data_vec <- as.numeric(unlist(truth_vector, use.names = FALSE))

fit_geometric <- function(data) {
  # 1) Coerce this chunk to a clean integer vector
  y <- as.numeric(unlist(data, use.names = FALSE))
  y <- y[!is.na(y)]               # drop NA
  y <- y[y > 0]                   # 1-based support: keep 1,2,...

  if (length(y) == 0) return(NA_real_)

  # 2) Ensure integer counts for geometric
  y <- round(y)
  x <- y - 1L                     # shift to 0-based for MASS::fitdistr("geometric")

  # 3) Fit 0-based geometric; 'prob' is the same p for 1-based geometric
  fit <- fitdistr(x, "geometric")
  unname(fit$estimate[["prob"]])
}

fit_geometric_sd <- function(data) {
  y <- as.numeric(unlist(data, use.names = FALSE))
  y <- y[!is.na(y)]
  y <- y[y > 0]

  if (length(y) == 0) return(NA_real_)

  y <- round(y)
  x <- y - 1L

  fit <- fitdistr(x, "geometric")
  # Standard error for 'prob' (guard against NULL)
  sd_prob <- fit$sd
  if (is.null(sd_prob)) return(NA_real_)
  unname(sd_prob[["prob"]])
}

# Split into 472 categories x 15 rows (adjust if your data length differs)
category_counts <- matrix(data_vec, nrow = 15, byrow = TRUE)

estimated_p_values <- apply(category_counts, 2, fit_geometric)
estimated_p_ses    <- apply(category_counts, 2, fit_geometric_sd)

# Drop NAs just in case (but KEEP endpoints 0 and 1)
estimated_p_values <- estimated_p_values[!is.na(estimated_p_values)]

## 2. Sample p per Species (hyperpriors) -----------------------

beta_a <- 1
beta_b <- 9

w_point <- 0.65
w_norm  <- 0.35
mu_norm <- 0.5
sd_norm <- 0.1

# helper: truncated normal on (0,1)
rtruncnorm <- function(n, mean, sd) {
  x <- rnorm(n, mean, sd)
  while (any(bad <- (x <= 0 | x >= 1))) {
    x[bad] <- rnorm(sum(bad), mean, sd)
  }
  x
}

# Unique Species = clusters
species_df <- gt %>% distinct(Species)

n_species <- nrow(species_df)

# Uninformed Beta(1,9)
p_beta <- rbeta(n_species, beta_a, beta_b)

# Mixture for informed prior
is_point <- rbinom(n_species, 1, w_point)
p_inf <- numeric(n_species)

p_inf[is_point == 1] <- 0.95

n_cont <- sum(is_point == 0)
if (n_cont > 0) {
  p_inf[is_point == 0] <- rtruncnorm(n_cont, mean = mu_norm, sd = sd_norm)
}

# Attach p-values back to each Species
species_ps <- species_df %>%
  mutate(
    p_beta = p_beta,
    p_inf  = p_inf
  )

## 3. Compute predicted k for EACH observation -----------------

gt_sim <- gt %>%
  left_join(species_ps, by = "Species") %>%
  mutate(
    k_obs  = Freq,
    k_beta = 1 + rgeom(n(), prob = p_beta),
    k_inf  = 1 + rgeom(n(), prob = p_inf)
  )

cat("Sanity check for k_obs, should show k=1 present:\n")
print(table(gt_sim$k_obs)[1:10])

## 4. Build long dataframe for k-plot --------------------------

long_df <- gt_sim %>%
  dplyr::select(k_obs, k_beta, k_inf) %>%  # avoid MASS::select clash
  pivot_longer(
    cols      = everything(),
    names_to  = "source",
    values_to = "k"
  ) %>%
  mutate(
    source = dplyr::recode(
      source,
      k_obs  = "Observed data",
      k_beta = "Prior predictive: Beta(1,9)",
      k_inf  = "Prior predictive: informed mixture"
    ),
    k = as.integer(k)
  ) %>%
  filter(k >= 1, k <= 20)

## 5. Long dataframe for p-plot (hyperpriors + data estimates) -

p_long_df <- bind_rows(
  tibble(
    p = p_beta,
    source = "Hyperprior: Beta(1,9)"
  ),
  tibble(
    p = p_inf,
    source = "Hyperprior: informed mixture"
  ),
  tibble(
    p = estimated_p_values,
    source = "Estimated from observed data"
  )
) %>%
  # KEEP p == 0 and p == 1; only drop obvious garbage
  filter(!is.na(p), p >= 0, p <= 1)

## 6. Shared colour mapping ------------------------------------

colour_map <- c(
  "Observed data"                      = "#1b9e77",
  "Prior predictive: uninformed (Beta(1,9))"        = "#7570b3",
  "Prior predictive: informed mixture" = "#d95f02",
  "Hyperprior:  uninformed (Beta(1,9)"              = "#7570b3",
  "Hyperprior: informed mixture"       = "#d95f02",
  "Estimated from observed data"       = "#1b9e77"  # match Observed data
)

## 7. PLOT: p-values panel (top) -------------------------------

p_pvals <- ggplot(p_long_df, aes(x = p)) +
  geom_histogram(
    aes(y = after_stat(density), fill = source),
    binwidth = 0.05,   # 20 bins: [0,0.05],...,[0.95,1.00]
    boundary = 0,      # align bins to start at 0
    colour   = "black",
    alpha    = 0.5
  ) +
  facet_wrap(~ source, ncol = 3) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.1)
  ) +
  coord_cartesian(xlim = c(0, 1)) +   # does NOT drop data before binning
  scale_fill_manual(values = colour_map) +
  labs(
    title = "Hyperpriors and data-based estimates of p (Geometric)",
    x     = "p",
    y     = "Density",
    fill  = NULL
  ) +
  theme_bw() +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90")
  )

## 8. PLOT: k-values panel (bottom) ----------------------------

p_kvals <- ggplot(long_df, aes(x = k)) +
  geom_histogram(
    aes(y = after_stat(density), fill = source),
    binwidth = 1,
    boundary = 0.5,     # bins centered on integers: [0.5,1.5], [1.5,2.5], ...
    colour   = "black",
    alpha    = 0.5
  ) +
  facet_wrap(~ source, ncol = 3) +
  scale_x_continuous(
    breaks = 1:20,
    limits = c(0.5, 20.5)   # ensure k=1 is fully inside first bin
  ) +
  scale_fill_manual(values = colour_map) +
  labs(
    title = "Observed vs prior predictive Geometric counts",
    x     = "k (Geometric count)",
    y     = "Density",
    fill  = NULL
  ) +
  theme_bw() +
  theme(
    legend.position  = "none",
    strip.background = element_rect(fill = "grey90")
  )

## 9. Combine panels and save ----------------------------------

combined_plot <- p_pvals / p_kvals  # p-panel on top, k-panel below

dir.create(file.path(repo_root, "Results/plots/priors"), recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path(repo_root, "Results/plots/priors/geom_priors_vs_observed_with_p_panel.jpeg"),
  plot   = combined_plot,
  width  = 10,
  height = 8,
  dpi    = 300
)
