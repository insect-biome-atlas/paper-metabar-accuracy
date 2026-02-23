###########################################################################
# R Plots for logZ
###########################################################################
library(readr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(jsonlite)
library(plotrix)

library(here)
repo_root <- here::here()

dir <- file.path(repo_root, "Results/step1_model_selection/smc_logZ/")

##############
optmodel15_comb_2theta <- fromJSON(paste0(dir,"combined_logZ/optmodel15_combined_2theta.json"))
############## Read in data for logZ plots
# Art
optmodel_art <- fromJSON(paste0(dir,"art_logZ/optmodel_art.json"))
optmodel15_1c <- fromJSON(paste0(dir,"art_logZ/optmodel15_1c.json"))
optmodel15_2k <- fromJSON(paste0(dir,"art_logZ/optmodel15_2k.json"))
optmodel15_2theta <- fromJSON(paste0(dir,"art_logZ/optmodel15_2theta.json"))
optmodel15_2k2theta <- fromJSON(paste0(dir,"art_logZ/optmodel15_2k2theta.json"))
# Bio
optmodel15_1c_biospikeins <- fromJSON(paste0(dir,"bio_logZ/optmodel15_1c_biospikeins.json"))
optmodel15_6k <- fromJSON(paste0(dir,"bio_logZ/optmodel15_6k.json"))
optmodel15_6theta <- fromJSON(paste0(dir,"bio_logZ/optmodel15_6theta.json"))
optmodel15_6k6theta <- fromJSON(paste0(dir,"bio_logZ/optmodel15_6k6theta.json"))
optmodel15_biospikeins <- fromJSON(paste0(dir,"bio_logZ/optmodel15_biospikeins.json"))
# Combined
optmodel15_8k <- fromJSON(paste0(dir,"combined_logZ/optmodel15_8k.json"))
optmodel15_8k8theta <- fromJSON(paste0(dir,"combined_logZ/optmodel15_8k8theta.json"))
optmodel15_1c_combined <- fromJSON(paste0(dir,"combined_logZ/optmodel15_1c_combined.json"))
optmodel15_combined <- fromJSON(paste0(dir,"combined_logZ/optmodel15_combined.json"))
# mcmc adapted models: fixed theta
optmodel15_smc_comb_fixTheta <- fromJSON(paste0(dir,"fix_theta_logZ/optmodel_smc_comb_fixTheta.json"))
optmodel15_smc_biospikeins_fixTheta <- fromJSON(paste0(dir,"fix_theta_logZ/optmodel_smc_biospikeins_fixTheta.json"))
optmodel15_art_fixTheta <- fromJSON(paste0(dir,"fix_theta_logZ/optmodel_smc_art_fixTheta.json"))
optmodel15_6k_fixTheta <- fromJSON(paste0(dir,"fix_theta_logZ/optmodel_smc_6k_fixTheta.json"))
optmodel15_8k_fixTheta <- fromJSON(paste0(dir,"fix_theta_logZ/optmodel_smc_8k_fixTheta.json"))

########################################################
#### Plot logZ
########################################################
jpeg(file.path(repo_root, "Results/plots/logZ/15samples - logZ.jpg"),
     width = 5,
     height = 9,
     units = "in",
     res = 300) 
par(mfrow = c(3, 1))
boxplot(
  optmodel15_1c$normConst,
  optmodel_art$normConst,
  optmodel15_art_fixTheta$normConst,
  optmodel15_2k$normConst,
  optmodel15_2theta$normConst,
  optmodel15_2k2theta$normConst,
  main = "A) Two synthetics spikeins",
  names = c(
    "1c",
    expression(paste("1k 1", theta)),
    expression(paste("1k ",theta,"=1")),
    expression(paste("2k 1", theta)),
    expression(paste("1k 2", theta)),
    expression(paste("2k 2", theta))
  ),
  col = adjustcolor(brewer.pal(6, "Blues"), alpha.f = 0.7),
  ylab = "log Z",
  outline=F
  #ylim = c(-680,-656) #c(-669,-656)
)
boxplot(
  optmodel15_1c_biospikeins$normConst,
  optmodel15_biospikeins$normConst,
  optmodel15_smc_biospikeins_fixTheta$normConst,
  optmodel15_6k$normConst,
  optmodel15_6k_fixTheta$normConst,
  optmodel15_6theta$normConst,
  optmodel15_6k6theta$normConst,
  main = "B) Six biological spikeins",
  names = c(
    "1c",
    expression(paste("1k 1", theta)),
    expression(paste("1k ",theta,"=1")),
    expression(paste("6k 1", theta)),
    expression(paste("6k ",theta,"=1")),
    expression(paste("1k 6", theta)),
    expression(paste("6k 6", theta))
  ),
  col = adjustcolor(brewer.pal(7, "Greens"), alpha.f = 0.7),
  ylab = "log Z",
  outline=F
  #ylim = c(-680,-656) #c(-669,-656)
)
boxplot(
  optmodel15_1c_combined$normConst,
  optmodel15_combined$normConst,
  optmodel15_smc_comb_fixTheta$normConst,
  optmodel15_comb_2theta$normConst,
  optmodel15_8k$normConst,
  optmodel15_8k_fixTheta$normConst,
  optmodel15_8k8theta$normConst,
  main = "C) Combined: two synthetic and six biological spikeins",
  names = c(
    "1c",
    expression(paste("3k 1", theta)),
    expression(paste("3k ",theta,"=1")),
    expression(paste("3k 2", theta)),
    expression(paste("8k 1", theta)),
    expression(paste("8k ",theta,"=1")),
    expression(paste("8k 8", theta))
  ),
  col = adjustcolor(brewer.pal(7, "Reds"), alpha.f = 0.7),
  ylab = "log Z",
  outline=F
  #ylim = c(-680,-656) #c(-669,-656)
)
dev.off()
##########################################  