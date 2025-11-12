# Compute the sd of relative error for lysate and homogenate
# Swedish IBA data, using each spikein as a calibration target,
# and 0 to 5 of the other spikeins for calibration

old_dir <- "~/Desktop/git/paper-metabar-accuracy/"
setwd("~/Desktop/git/utils_Nov2025/")
source("get_iba_co1_data_fxn.R")
setwd(old_dir)


# Define a function that generates all possible bool vectors
generate_includes <- function(len) {

    includes <- rep(FALSE,times=len)
    include_list <- list(includes)
    
    includes <- increment(includes)
    while (sum(includes)!=0) {
        include_list[[length(include_list)+1]] <- includes
        includes <- increment(includes)
    }

    return (include_list)
}

# Define a function that increments a bool vector
increment <- function(includes) {

    for (i in length(includes):1) {
        if (includes[i]==FALSE)
            break
    }
    if (includes[i]==TRUE)
        j <- i
    else {
        includes[i] <- TRUE
        j <- i+1
    }
    while (j<=length(includes)) {
        includes[j] <- FALSE
        j <- j + 1
    }
    return (includes)
}


# Get homogenate and lysate data
data_path <- "~/Desktop/git/figshare-repos/iba/processed_data/v3/"
metadata_path <- "~/Desktop/git/figshare-repos/iba/raw_data/v6/"
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

# Make homogenate data version with only bio spikeins
HSAB <- HSA[!grepl("synth",HSA$cluster),]

# List names of spikeins
bio_spikes <- c(
                "Blattidae_cluster1",
                "Gryllidae_cluster1",
                "Gryllidae_cluster2",
                "Drosophilidae_cluster1",
                "Drosophilidae_cluster2",
                "Drosophilidae_cluster3"
                )

# Define specimen numbers for spikeins
bio_specimens <- c(2,1,1,3,1,1)

# Define x axis labels
x_labels <- c("Sh.la","Gr.bi","Gr.su","Dr.bi","Dr.se","Dr.ja")

# Define function for computing errors (deviations) with or without
# calibration with various number of calibrator species
compute_cal_error <- function(D) {

    # Set up result df
    res <- data.frame(list(
                    spikein=character(),
                    num_includes=numeric(),
                    Sh.la=logical(),
                    Gr.bi=logical(),
                    Gr.su=logical(),
                    Dr.bi=logical(),
                    Dr.se=logical(),
                    Dr.ja=logical(),
                    stdev_log10=numeric(),
                    stdev_rel=numeric(),
                    iqr_log10=numeric(),
                    iqr_rel=numeric()
                    )
                )

    # Iterate over biological spikeins
    for (i in 1:length(bio_spikes)) {

        spikein <- bio_spikes[i]
        x <- as.numeric(D[i,-1] / bio_specimens[i])  # Get the target reads per specimen

        # Get all possible bool vectors
        include_list <- generate_includes(6)

        for (j in 1:length(include_list)) {

            includes <- include_list[[j]]
            
            # Skip this bool vector if the target is included
            if (includes[i]==TRUE)
                next

            num_includes <- sum(includes)
            
            # Use 0 includes to get uncalibrated coefficient of variation, otherwise
            # we compute the coefficient of variation for calibrated reads
            if (num_includes==0) {
                stdev_log10 <- sd(log10(x/mean(x)))
                iqr_log10 <- IQR(log10(x/mean(x)))
            } else {
                y <- as.numeric(colSums(D[includes,-1]))
                fit <- lm(x ~ y + 0)
                stdev_log10 <- sd(log10(x/fit$fitted.values))
                iqr_log10 <- IQR(log10(x/fit$fitted.values))
            }
            res <- rbind(res,
                   list(spikein=spikein,
                        num_includes=num_includes,
                        Sh.la=includes[1],
                        Gr.bi=includes[2],
                        Gr.su=includes[3],
                        Dr.bi=includes[4],
                        Dr.se=includes[5],
                        Dr.ja=includes[6],
                        stdev_log10=stdev_log10,
                        stdev_rel=10^stdev_log10,
                        iqr_log10=iqr_log10,
                        iqr_rel=10^iqr_log10
                        )
                   )
        }
    }

    return(res)
}

# Compute lysate and homogenate coefficients of variation
cat("Computing relative errors for lysate data\n")
lysate_cal_error <- compute_cal_error(LSA)
cat("Computing relative errors for homogenate data\n")
homogenate_cal_error <- compute_cal_error(HSAB)

