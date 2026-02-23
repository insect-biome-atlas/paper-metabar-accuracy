# run_postprocessing_after_mcmc.R

cat("Starting post-processing after MCMC...\n")

library(jsonlite)
library(here)
repo_root <- here::here()
source(file.path(repo_root, "R/pre-post-processing/postprocessing_after_mcmc.R"))

dir_in     <- file.path(repo_root, "Results/step1_model_selection/mcmc/chains/")
dir_out    <- file.path(repo_root, "Results/step2_predictions/prediction_inputs/spikeins/")
prior_data <- fromJSON(file.path(repo_root, "Data/15_hom_priors.json"), simplifyVector = FALSE)

############## Read in data from mcmc
post_data_art <- fromJSON(paste0(dir_in,"optmodel_mcmc_art_fixTheta_1.json"))
post_data_bio <- fromJSON(paste0(dir_in,"optmodel_mcmc_biospikeins_fixTheta_1.json"))
post_data_comb <- fromJSON(paste0(dir_in,"optmodel_mcmc_comb_fixTheta_1.json"))
post_data_6k <- fromJSON(paste0(dir_in,"optmodel_mcmc_6k_fixTheta_1.json"))
post_data_8k <- fromJSON(paste0(dir_in,"optmodel_mcmc_8k_fixTheta_1.json"))

############## Post-processing step 
res_art <- post_processing_single_mcmc(post_data_art)
res_bio <- post_processing_single_mcmc(post_data_bio)
res_comb <- post_processing_comb_single_mcmc(post_data_comb)
res_6k <- post_processing_comb_6k_mcmc(post_data_6k)
res_8k <- post_processing_comb_8k_mcmc(post_data_8k)

############## Output for prediction
output_for_pred_single_mcmc(res_art, prior_data, paste0(dir_out,"pred_in_art_mcmc.json"))
output_for_pred_single_mcmc(res_bio, prior_data, paste0(dir_out,"pred_in_bio_mcmc.json"))
output_for_pred_single_mcmc(res_comb, prior_data, paste0(dir_out,"pred_in_comb_mcmc.json"))
output_for_pred_comb_mcmc(res_6k, prior_data, paste0(dir_out,"pred_in_6k_mcmc.json"))
output_for_pred_comb_mcmc(res_8k, prior_data, paste0(dir_out,"pred_in_8k_mcmc.json"))

############## Finished
cat("Done!\n")