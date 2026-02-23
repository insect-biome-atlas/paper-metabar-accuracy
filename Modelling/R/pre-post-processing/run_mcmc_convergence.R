#Usage: source("treeppl_mcmc_convergence.R")
# Below for the models with fixed theta

library(here)
repo_root  <- here::here()
chains_dir <- file.path(repo_root, "Results/step1_model_selection/mcmc/chains")
conv_dir   <- file.path(repo_root, "Results/step1_model_selection/mcmc/convergence")
source(file.path(repo_root, "R/pre-post-processing/mcmc_convergence.R"))
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

run_convergence(
  chain_paths = file.path(chains_dir, c(
    "optmodel_mcmc_biospikeins_fixTheta_1.json",
    "optmodel_mcmc_biospikeins_fixTheta_2.json",
    "optmodel_mcmc_biospikeins_fixTheta_3.json"
  )),
  burn_in = 0.5,            # proportion or absolute integer; 0.5 means drop first 50%
  thin = 1,                 # keep every 'thin' sample
  rhat_thresh = 1.05,       # threshold for convergence decision
  ess_thresh = 300,         # effective sample size threshold
  out_dir = file.path(conv_dir, "mcmc_check_biospikeins")
)

run_convergence(
  chain_paths = file.path(chains_dir, c(
    "optmodel_mcmc_comb_fixTheta_1.json",
    "optmodel_mcmc_comb_fixTheta_2.json",
    "optmodel_mcmc_comb_fixTheta_3.json"
  )),
  burn_in = 0.5,            # proportion or absolute integer; 0.5 means drop first 50%
  thin = 1,                 # keep every 'thin' sample
  rhat_thresh = 1.05,       # threshold for convergence decision
  ess_thresh = 300,         # effective sample size threshold
  out_dir = file.path(conv_dir, "mcmc_check_comb")
)

run_convergence(
  chain_paths = file.path(chains_dir, c(
    "optmodel_mcmc_6k_fixTheta_1.json",
    "optmodel_mcmc_6k_fixTheta_2.json",
    "optmodel_mcmc_6k_fixTheta_3.json"
  )),
  burn_in = 0.5,            # proportion or absolute integer; 0.5 means drop first 50%
  thin = 1,                 # keep every 'thin' sample
  rhat_thresh = 1.05,       # threshold for convergence decision
  ess_thresh = 300,         # effective sample size threshold
  out_dir = file.path(conv_dir, "mcmc_check_6k")
)

run_convergence(
  chain_paths = file.path(chains_dir, c(
    "optmodel_mcmc_8k_fixTheta_1.json",
    "optmodel_mcmc_8k_fixTheta_2.json",
    "optmodel_mcmc_8k_fixTheta_3.json"
  )),
  burn_in = 0.5,            # proportion or absolute integer; 0.5 means drop first 50%
  thin = 1,                 # keep every 'thin' sample
  rhat_thresh = 1.05,       # threshold for convergence decision
  ess_thresh = 300,         # effective sample size threshold
  out_dir = file.path(conv_dir, "mcmc_check_8k")
)

run_convergence(
  chain_paths = file.path(chains_dir, c(
    "optmodel_mcmc_art_fixTheta_1.json",
    "optmodel_mcmc_art_fixTheta_2.json",
    "optmodel_mcmc_art_fixTheta_3.json"
  )),
  burn_in = 0.3,            # proportion or absolute integer; 0.5 means drop first 50%
  thin = 1,                 # keep every 'thin' sample
  rhat_thresh = 1.05,       # threshold for convergence decision
  ess_thresh = 300,         # effective sample size threshold
  out_dir = file.path(conv_dir, "mcmc_check_art")
)
