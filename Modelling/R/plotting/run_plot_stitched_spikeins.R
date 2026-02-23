
library(here)
repo_root <- here::here()
source(file.path(repo_root, "R/plotting/plot_from_stitched_spikeins.R"))

spike_names_vec <- c(
  "Shelfordella lateralis (2)",
  "Drosophila bicornuta (3)",
  "Drosophila serrata (1)",
  "Drosophila jambulina (1)",
  "Gryllus bimaculatus (1)",
  "Gryllodes sigillatus (1)"
)
true_n_vec      <- c(2,3,1,1,1,1)

#base_dir <- file.path(repo_root, "Results/step2_predictions/stitched/spikeins/uninform")
base_dir <- file.path(repo_root, "Results/step2_predictions/stitched/spikeins/inform_pfit")
#base_dir <- file.path(repo_root, "Results/step2_predictions/stitched/spikeins/inform_adjust")

models <- list(
  combined = file.path(base_dir, "stitched_pred_spikeins_comb.rds"),
  bio      = file.path(base_dir, "stitched_pred_spikeins_bio.rds"),
  #art      = file.path(base_dir, "stitched_pred_spikeins_art.rds"),
  match6k  = file.path(base_dir, "stitched_pred_spikeins_6k.rds"),
  match8k  = file.path(base_dir, "stitched_pred_spikeins_8k.rds")
)

# check?
#print(models)
#file.exists(unlist(models))

run_plotter_multi(
  models = models,
  output_dir = file.path(repo_root, "Results/plots/spikeins/inform_pfit/"),
  spikein_names_src = paste(spike_names_vec, collapse = ","),  
  true_n_src        = paste(true_n_vec, collapse = ","),        
  x_max = 12,
  agg_mode       = "all"                  # all | mean | map | wmean_normconst
)
