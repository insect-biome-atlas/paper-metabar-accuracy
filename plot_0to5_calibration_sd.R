# Plot the log read sd for lysate and homogenate
# Swedish IBA data and 0 to 5 spike-in calibrations

# Libraries needed
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

old_dir <- getwd()
setwd("~/dev/ms-repos-iba/utils/")
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
cat("Getting homogenate data\n"
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

cat("Computing plots\n")


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

# Define x axis labels
x_labels <- c("Sh.la","Gr.bi","Gr.su","Dr.bi","Dr.se","Dr.ja")


# Compute calibration sd for L counts with 0 to 5 spike-ins
# -------------------------------------------------------------------

# Set up data frames
D <- data.frame(list(
                    spikein=character(),
                    num_includes=numeric(),
                    Sh.la=logical(),
                    Gr.bi=logical(),
                    Gr.su=logical(),
                    Dr.bi=logical(),
                    Dr.se=logical(),
                    Dr.ja=logical(),
                    coeff_var=numeric()
                    )
                )
E <- D

# Iterate over biological spikeins
for (i in 1:length(bio_spikes)) {

    spikein <- bio_spikes[i]
    x <- log10(LSA[i,-1])

    # Get all possible bool vectors
    include_list <- generate_includes(6)

    for (j in 1:length(include_list)) {

        includes <- include_list[[j]]
        if (includes[i]==TRUE)
            next

        num_includes <- sum(includes)
        if (num_includes==0) {
            stddev <- sd(x)
        } else {
            y <- log10(colSums(LSA[includes,-1]))
            f <- y - mean(y)
            stddev <- sd(x - f)
        }
        D <- rbind(D,
                   list(spikein=spikein,
                        num_includes=num_includes,
                        Sh.la=includes[1],
                        Gr.bi=includes[2],
                        Gr.su=includes[3],
                        Dr.bi=includes[4],
                        Dr.se=includes[5],
                        Dr.ja=includes[6],
                        stddev=stddev
                        )
                   )
    }
}


# Compute calibration sd for H counts with 0 to 5 spike-ins
# -------------------------------------------------------------------

# Iterate over biological spikeins
for (i in 1:length(bio_spikes)) {

    spikein <- bio_spikes[i]
    x <- log10(HSAB[i,-1])

    # Get all possible bool vectors
    include_list <- generate_includes(6)

    for (j in 1:length(include_list)) {

        includes <- include_list[[j]]
        if (includes[i]==TRUE)
            next
 
        num_includes <- sum(includes)
        if (num_includes==0) {
            stddev <- sd(x)
        } else {
            y <- log10(colSums(HSAB[includes,-1]))
            f <- y - mean(y)
            stddev <- sd(x - f)
        }
        E <- rbind(E,
                   list(spikein=spikein,
                        num_includes=num_includes,
                        Sh.la=includes[1],
                        Gr.bi=includes[2],
                        Gr.su=includes[3],
                        Dr.bi=includes[4],
                        Dr.se=includes[5],
                        Dr.ja=includes[6],
                        stddev=stddev
                        )
                   )
    }
}


# Make plots
p1 <- ggplot(data=D, aes(x=num_includes, y=coeff_var, colour=factor(spikein)) +
             theme_minimal() +
             scale_x_discrete(name="No. calibrator spike-ins", labels=x_labels) +
             geom_point(size=3) +
             scale_y_continuous(name="Coefficient of variation",trans="log10")
         )
p2 <- ggplot(data=E, aes(x=num_includes, y=coeff_var, colour=factor(spikein)) +
             theme_minimal() +
             scale_x_discrete(name="No. calibrator spike-ins", labels=x_labels) +
             geom_point(size=3) +
             scale_y_continuous(name="Coefficient of variation",trans="log10")
         )
ggsave(file="Figures/Fig_0to5_calibrations_CV.jpg", height=3.5, width=7, plot = (p1 + p2) + plot_layout(axis_titles="collect",ncol=2) + plot_annotation(tag_levels="A"))


#––– Single-spikein performance –––––––––––––––––––––––––––––––––––––––––––––––––

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

p_simple <- ggplot(plot_df,
                   aes(x = spikein, y = stddev,
                       shape = calibrator, colour = calibrator)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_y") +
  scale_shape_discrete(drop = FALSE) +
  scale_colour_manual(values = custom_cols, drop = FALSE) +
  labs(
    title  = "Residual SD (log10) for each target spikein\n(single-spike-in calibrations)",
    x      = "Target spike-in                                                         Target spike-in",
    y      = "Residual stddev",
    shape  = "Calibrator",
    colour = "Calibrator"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("Fig_0to5_single_spikein_calibration.pdf",
       plot = p_simple, width = 8, height = 4)

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
                           colour = calibrator, shape = calibrator)) +
  geom_jitter(width = 0.002, height = 0, size = 3, alpha = 0.8) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = custom_cols, drop = FALSE) +
  scale_shape_discrete(drop = FALSE) +
  labs(
    title  = "Calibration efficiency\nresidual SD of target vs calibrator read-proportion",
    x      = "Proportion of total reads of the calibrator spike-in",
    y      = "Residual SD (log10) of the target spike-in",
    colour = "Calibrator",
    shape  = "Calibrator"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave("Fig_0to5_spikein_efficiency_by_calibrator.pdf",
       plot = p_eff_by_cal, width = 8, height = 4)
