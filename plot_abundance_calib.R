# Plot the lognormal read variance for various
# Swedish IBA data and spike-in calibrations

library(reshape2)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(patchwork)

old_dir <- getwd()
setwd("~/dev/ms-repos-iba/utils/")
source("get_iba_co1_data_fxn.R")
setwd(old_dir)

# List names of spikeins
bio_spikes <- c(
                "Blattidae_cluster1",
                "Gryllidae_cluster1",
                "Gryllidae_cluster2",
                "Drosophilidae_cluster1",
                "Drosophilidae_cluster2",
                "Drosophilidae_cluster3"
                )
synth_spikes <- c("Callio-synth","tp53-synth")
spike_species <- c("Shelfordella lateralis","Gryllus bimaculatus","Gryllodes sigillatus ","Drosophila bicornuta","Drosophila serrata","Drosophila jambolina")

# Define specimen numbers for spikeins
bio_specimens <- c(2,1,1,3,1,1)

# Define x axis labels
x_labels <- c("Sh.la","Gr.bi","Gr.su","Dr.bi","Dr.se","Dr.ja")

# Define plot function 
vio_plot <- function(D, x_labels, plot_title = NULL) {
  ggplot(data = D, aes(x = cluster, y = log10(reads), fill = cluster)) +
    theme_minimal() +
    geom_violin(alpha = 0.6, color = "grey30") +
    ylim(c(-0.05, 6.5)) +
    scale_x_discrete(name = NULL, labels = NULL) +  # drop x labels
    scale_fill_manual(values = custom_cols, labels = spike_species, name = "Species") +
    ylab("Reads per specimen (log10)") +
    ggtitle(plot_title) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

# Define calibration functions
calibrate_reads_bio <- function(X) {

    # Convert input data to long format
    D <- reshape2::melt(X,id="cluster",variable.name="sample",value.name="reads")

    # Leave-One-Out (LOO) calibration
    X_LOO <- X[FALSE,]
    for (i in 1:nrow(X)) {
        x <- colSums(X[-i,-1])
        X_LOO <- rbind(X_LOO,as.list(c(cluster=X$cluster[i],x)))
    }

    E <- reshape2::melt(X_LOO,id="cluster",variable.name="sample",value.name="spikein_reads")

    D$cluster <- factor(D$cluster,levels=bio_spikes)
    E$cluster <- factor(E$cluster,levels=bio_spikes)
    D$reads <- as.numeric(D$reads)
    E$spikein_reads <- as.numeric(E$spikein_reads)
    D <- merge(D,E)

    D$mean <- numeric(nrow(D))
    D$residual <- numeric(nrow(D))

    # Fit the linear model
    for (i in 1:length(bio_spikes)) {
        cluster <- bio_spikes[i]
        specimens <- bio_specimens[i]
        idx <- D$cluster == cluster
        D$reads[idx] <- D$reads[idx] / specimens    # Measure reads/specimen
        fit <- lm(D$reads[idx]~D$spikein_reads[idx]+0)
        D$residual[idx] <- fit$residuals
        D$fitted.values[idx] <- fit$fitted.values
        D$mean[idx] <- mean(D$reads[idx])
    }
    D$calibrated_reads <- mean(fit$fitted.values)*(D$reads/fit$fitted.values)
    
    return(D)
}

calibrate_reads_bio_synth <- function(X) {

    # Make data version with only synthetic spikeins
    XS <- X[grepl("synth",X$cluster),]
    
    # Make data version with only biological spikeins
    XB <- X[!grepl("synth",X$cluster),]
    
    # Calibrate with only bio spikeins
    D <- calibrate_reads_bio(XB)

    # Rename spikein_reads to bio_reads
    idx <- which(colnames(D)=="spikein_reads")
    colnames(D)[idx] <- "bio_reads"

    # Rename calibrated_reads to bio_calibrated_reads
    idx <- which(colnames(D)=="calibrated_reads")
    colnames(D)[idx] <- "bio_calibrated_reads"

    # Summarize synthetic spike-in reads
    XB_S <- XB[FALSE,]
    x <- colSums(XS[,-1])
    for (i in 1:nrow(XB)) {
        XB_S <- rbind(XB_S,as.list(c(cluster=XB$cluster[i],x)))
    }
    E <- reshape2::melt(XB_S,id="cluster",variable.name="sample",value.name="synth_reads")

    E$cluster <- factor(E$cluster,levels=bio_spikes)
    E$synth_reads <- as.numeric(E$synth_reads)

    D <- merge(D,E)

    # Fit the linear models for synth and bio+synth spikeins
    for (i in 1:length(bio_spikes)) {
        cluster <- bio_spikes[i]
        specimens <- bio_specimens[i]
        idx <- D$cluster == cluster
        fit1 <- lm(D$reads[idx]~D$synth_reads[idx]+0)
        fit2 <- lm(D$reads[idx]~D$synth_reads[idx]+D$bio_reads[idx]+0)
#        print(summary(fit2))
        D$residual1[idx] <- fit1$residuals
        D$residual2[idx] <- fit2$residuals
        D$fitted.values1[idx] <- fit1$fitted.values
        D$fitted.values2[idx] <- fit2$fitted.values
    }
    D$synth_calibrated_reads <- mean(fit1$fitted.values)*(D$reads/fit1$fitted.values)
    D$spikein_calibrated_reads <- mean(fit2$fitted.values)*(D$reads/fit2$fitted.values)

    return(D)
}

# Get homogenate and lysate data
data_path <- "~/dev/figshare-repos/iba/processed_data/v3/"
metadata_path <- "~/dev/figshare-repos/iba/raw_data/v6/"
cat("Getting lysate data\n")
L <- get_iba_co1_data(
                      data_path=data_path,
                      metadata_path=metadata_path,
                      country="SE",
                      dataset="lysate",
                      calibrate=FALSE,
                      remove_spikes=FALSE
                    )
cat("Getting homogenate data\n")
H <- get_iba_co1_data(
                      data_path=data_path,
                      metadata_path=metadata_path,
                      country="SE",
                      dataset="homogenate",
                      calibrate=FALSE,
                      remove_spikes=FALSE
                     )
cat("Computing plots\n")

L <- data.frame(L)
H <- data.frame(H)

# Find sample columns
colnames(L) -> x
sample_cols_L <- which(grepl("P",x) & grepl("_",x))
colnames(H) -> x
sample_cols_H <- which(grepl("P",x) & grepl("_",x))

# Find spikeins
prop_samples_L <- rowMeans(L[,sample_cols_L]>0)
prop_samples_H <- rowMeans(H[,sample_cols_H]>0)

spikeins_L <- which(prop_samples_L > 0.95)
spikeins_H <- which(prop_samples_H > 0.95)

LS <- L[spikeins_L,c(which(colnames(L)=="cluster"),sample_cols_L)]
HS <- H[spikeins_H,c(which(colnames(H)=="cluster"),sample_cols_H)]

# Only keep columns containing all spikeins
LSA <- LS[,c(TRUE,colMeans(LS[,-1]>0)>0.95)]
HSA <- HS[,c(TRUE,colMeans(HS[,-1]>0)>0.95)]


# Violin plots for L counts (raw/calibrated with biological spike-ins)
# -------------------------------------------------------------------

# Calibrate reads based on bio spikeins; this will return data in long format
D <- calibrate_reads_bio(LSA)

# Generate the plots
p1 <- vio_plot(D, x_labels, plot_title = "Uncalibrated")
p2 <- vio_plot(data.frame(list(cluster=D$cluster, reads=D$calibrated_reads)), x_labels, plot_title = "Bio spike-in calibrated")
# ggsave(file="Fig_lysate_cals.jpg", height=3.5, width=7, plot = p1 + p2)
p1
p2
# Violin plots for H counts (raw/calibrated using syn/bio/both spike-ins) 
# -----------------------------------------------------------------------

# Calibrate reads based on bio and synth spikeins; this will return data in long format
DH <- calibrate_reads_bio_synth(HSA)

p3 <- vio_plot(DH, x_labels, plot_title = "Uncalibrated")
p4 <- vio_plot(data.frame(list(cluster=DH$cluster, reads=DH$synth_calibrated_reads)), x_labels, plot_title = "Synth spike-in calibrated")
p5 <- vio_plot(data.frame(list(cluster=DH$cluster, reads=DH$bio_calibrated_reads)), x_labels, plot_title = "Bio spike-in calibrated")
p6 <- vio_plot(data.frame(list(cluster=DH$cluster, reads=DH$spikein_calibrated_reads)), x_labels, plot_title = "Combined spike-ins calibrated")
# ggsave(file="Fig_homogenate_cals.jpg", height=7, width=7, plot = (p3 + p4)/(p5 + p6))

combined_plot <- (p1 + p2 + p3 + p4 + p5 + p6) +
  plot_layout(ncol = 2, guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")  

ggsave(file = "Figures/Fig_calibrations.jpg",
       plot = combined_plot,  height = 10.5,  width = 7,  dpi = 300)

# Repeat plots with only samples for which there is both lysate and homogenate data
# ---------------------------------------------------------------------------------

# Restrict LSA and HSA datasets to shared samples
M <- read.delim("~/dev/figshare-repos/iba/raw_data/v6/CO1_sequencing_metadata_SE.tsv")
field_sample <- function(x) { M$sampleID_FIELD[match(x,M$sampleID_NGI)] }
x <- field_sample(colnames(LSA)[-1])
y <- field_sample(colnames(HSA)[-1])
LSA_idx <- c(1, 1+which(x %in% y))
HSA_idx <- c(1, 1+which(y %in% x))
LSS <- LSA[,LSA_idx]
HSS <- HSA[,HSA_idx]

# Calibrate the reads
D <- calibrate_reads_bio(LSA)
DH <- calibrate_reads_bio_synth(HSA)

# Generate the plots
p1 <- vio_plot(D, x_labels, plot_title = "Uncalibrated")
p2 <- vio_plot(data.frame(list(cluster=D$cluster, reads=D$calibrated_reads)), x_labels, plot_title = "Bio spike-in calibrated")
p3 <- vio_plot(DH, x_labels, plot_title = "Uncalibrated")
p4 <- vio_plot(data.frame(list(cluster=DH$cluster, reads=DH$synth_calibrated_reads)), x_labels, plot_title = "Synth spike-in calibrated")
p5 <- vio_plot(data.frame(list(cluster=DH$cluster, reads=DH$bio_calibrated_reads)), x_labels, plot_title = "Bio spike-in calibrated")
p6 <- vio_plot(data.frame(list(cluster=DH$cluster, reads=DH$spikein_calibrated_reads)), x_labels, plot_title = "Combined spike-ins calibrated")

combined_plot_shared <- (p1 + p2 + p3 + p4 + p5 + p6) +
  plot_layout(ncol = 2, guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = "A") & theme(legend.position = "bottom")  

# Join all calibrations for shared samples in one figure
ggsave(file = "Figures/Fig_calibrations_shared.jpg",
       plot = combined_plot_shared,  height = 10.5,  width = 7,  dpi = 300)

# Check if calibrated reads are approximately log_normal
# ------------------------------------------------------

# Set spikein_names
spikein_names <- c("Shelfordella lateralis",
                   "Gryllus bimaculatus",
                   "Gryllodes sigillatus",
                   "Drosophila bicornuta",
                   "Drosophila serrata",
                   "Drosopila jambolina")

# Spikein-calibrated lysate reads
dens_plot <- function (df, spikein, spikein_name) {
    df <- df[df$cluster==spikein,]
    df$log10_rel_error <- log10(df$reads/df$fitted.values)
    ggplot(df, aes(x = log10_rel_error)) +
        geom_density(fill = "steelblue", alpha = 0.4) +
        labs(title = spikein_name,
            x = "Relative error (log10)",
            y = "Density") +
        theme_minimal()
}

p1 <- dens_plot(D, bio_spikes[1], spikein_names[1])
p2 <- dens_plot(D, bio_spikes[2], spikein_names[2])
p3 <- dens_plot(D, bio_spikes[3], spikein_names[3])
p4 <- dens_plot(D, bio_spikes[4], spikein_names[4])
p5 <- dens_plot(D, bio_spikes[5], spikein_names[5])
p6 <- dens_plot(D, bio_spikes[6], spikein_names[6])

ggsave(file="Figures/Fig_rel_cal_error_lysates.jpg", height=10.5, width=7, plot = (p1 + p2 + p3 + p4 + p5 + p6) + plot_layout(axis_titles="collect",ncol=2) + plot_annotation(tag_levels="A"))

# Bio spike-in-calibrated homogenate reads
p1 <- dens_plot(DH, bio_spikes[1], spikein_names[1])
p2 <- dens_plot(DH, bio_spikes[2], spikein_names[2])
p3 <- dens_plot(DH, bio_spikes[3], spikein_names[3])
p4 <- dens_plot(DH, bio_spikes[4], spikein_names[4])
p5 <- dens_plot(DH, bio_spikes[5], spikein_names[5])
p6 <- dens_plot(DH, bio_spikes[6], spikein_names[6])

ggsave(file="Figures/Fig_rel_cal_error_homogenates.jpg", height=10.5, width=7, plot = (p1 + p2 + p3 + p4 + p5 + p6) + plot_layout(axis_titles="collect",ncol=2) + plot_annotation(tag_levels="A"))
