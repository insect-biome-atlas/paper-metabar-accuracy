# Plot matches of lysate and homogenate reads

library(reshape2)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(stringr)

# Get matched data
D <- read.delim("matched_lysate_homogenate.tsv")
D <- D[D$lysate_reads >0 | D$homogenate_reads > 0,]

#EXPLORATORY
## limiting number of taxa shown
top_taxa <- D %>%
  count(Order, sort = TRUE) %>%
  slice_head(n = 10) %>%
  pull(Order)
D_top10orders <- D %>% filter(Order %in% top_taxa)

ggplot(D_top10orders, aes(x = lysate_reads, y = homogenate_reads, color = Order)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Lysate Reads (log10)",
    y = "Homogenate Reads (log10)",
    title = "Scatter Plot of Lysate vs. Homogenate Reads (Log Scale)",
    color = "Taxonomy"
  ) +
  theme_minimal() +
  theme(legend.position = "right")


#Number of samples that were both lysed and homogenized?
length(unique(D$sample))

# Look at reads per sample for homogs and lysates:
# Sum reads per sample and treatment
D_read_summary <- D %>%
  group_by(sample) %>%
  summarize(
    lysate_total = sum(lysate_reads, na.rm = TRUE),
    homogenate_total = sum(homogenate_reads, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(lysate_total, homogenate_total),
               names_to = "treatment",
               values_to = "total_reads") %>%
  mutate(treatment = ifelse(treatment == "lysate_total", "Lysate", "Homogenate"))

## filter out anything but Arthropoda
D_filtered <- D %>%
  filter(Phylum %in% "Arthropoda")
length(unique(D_filtered$cluster))
length(unique(D_filtered$sample))
D_filtered$taxon <- paste0(D_filtered$Order, "_", D_filtered$Family)
D_filtered$taxon2 <- paste0(D_filtered$Family, " (", D_filtered$Order, ")")

D_family <- D_filtered %>%
  group_by(Order, Family) %>%   # keep both order + family
  summarise(
    lysate_mean = mean(lysate_reads, na.rm = TRUE),
    homogenate_mean = mean(homogenate_reads, na.rm = TRUE),
    lysate_sd = sd(lysate_reads, na.rm = TRUE),
    homogenate_sd = sd(homogenate_reads, na.rm = TRUE),
    .groups = "drop"
  )

plot<- ggplot(D_family, aes(x = lysate_mean + 1, y = homogenate_mean + 1, color = Order)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_errorbar(aes(ymin = homogenate_mean + 1 - homogenate_sd,
                    ymax = homogenate_mean + 1 + homogenate_sd),
                width = 0.1, alpha = 0.4) +
  geom_errorbarh(aes(xmin = lysate_mean + 1 - lysate_sd,
                     xmax = lysate_mean + 1 + lysate_sd),
                 height = 0.1, alpha = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Lysate Reads (+1)",
    y = "Homogenate Reads (+1)",
    title = "Family-level means ± SD (colored by Order)",
    color = "Order"
  ) +
  theme_minimal()
ggsave(file="Fig_S4_Scatter_all_famillies.jpg", height=7, width=14, plot = plot)
plot

########## Calculate Lysate-to-Homogenate RATIO
# average read depth varies between lysate and homogenate. find out by how much:
D_lys_sample_mean <- sum(D$lysate_reads)/length(unique(D$sample))
D_hom_sample_mean <- sum(D$homogenate_reads)/length(unique(D$sample))
factor<-D_lys_sample_mean/D_hom_sample_mean

D$homogenate_reads_corrected <- D$homogenate_reads * factor
#checking id the means of "lysate reads" and "homogenate_reads_corrected" are the same.
summary(D)
#they are. we proceed

## Exclude all pairs of samples that have 0 in any treatment.
#and calculate the ratio of lysate reads to (corrected) homogenate reads.
D$lys_homog_ratio <- ifelse(
  D$lysate_reads != 0 & D$homogenate_reads != 0,
  D$lysate_reads / D$homogenate_reads_corrected,
  NA  # Use NA for rows where one or both values are zero
)

#Filter anything that is not Arthropoda
D_filtered <- D %>%
  filter(Phylum %in% "Arthropoda")
length(unique(D_filtered$cluster))
D_filtered$taxon <- paste0(D_filtered$Order, "_", D_filtered$Family)
D_filtered$taxon2 <- paste0(D_filtered$Family, " (", D_filtered$Order, ")")


#Filter out the NAs for plotting
D_filtered_plot <- D_filtered[!is.na(D_filtered$lys_homog_ratio) & is.finite(D_filtered$lys_homog_ratio), ]
#Additionally, eliminate all taxa that are unclassifies (which is essentially few Diptera clusters)
D_filtered_plot <- D_filtered_plot %>%
  filter(!str_detect(taxon2, "unclassified"))
D_filtered_plot <- D_filtered_plot %>%
  filter(!str_detect(taxon2, "_X"))


# Filter taxa with at least 100 data points
D_filtered_plot_500 <- D_filtered_plot %>%
  group_by(taxon) %>%
  filter(n() >= 500) %>%
  ungroup()


# Boxplot split by Order
as<-ggplot(D_filtered_plot_500, aes(x = lys_homog_ratio, y = Family, fill = Order)) +
  geom_violin(alpha = 0.5, trim = TRUE, color = NA) +   # density shape
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black") +
  #geom_jitter(width = 0.05, size = 0.7, alpha = 0.3, color = "black") +  # small & transparent points
  #geom_jitter(
  # data = ~ dplyr::sample_n(., min(500, nrow(.))),  # max 500 points per facet
  #  width = 0.05, size = 0.7, alpha = 0.3
  #) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_log10() +
  facet_grid(rows = vars(Order), scales = "free_y", space = "free_y") +
  theme_minimal() +
  labs(
    x = "Lysate / Homogenate Read Ratio (log scale)",
    y = "",
    title = ""
  ) +
  theme(
    legend.position = "none",
    strip.text.y = element_text(size = 12, angle = 0), # Order labels, ie face = "bold",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "gray90", color = "black"),
    axis.text.y = element_text(size = 12)   # Family labels 
  )
as
ggsave(file="Fig2_All_Samples_500_Ratio.jpg", height=13, width=7, plot = as)

### "ZERO clusters"
####### Look at those clusters that had count of 0 in one or the other treatment (lysis/homog)
# subset with zeros in either lysate or homogenate
#D2 <- D_filtered %>%
# mutate(
#  lysate_reads = ifelse(lysate_reads < 20, 0, lysate_reads),
#  homogenate_reads = ifelse(homogenate_reads < 20, 0, homogenate_reads)
#)

# Get matched data again (cause RATIO plotting turned all rows with "0" into NAs)
D2 <- read.delim("matched_lysate_homogenate.tsv")
D2 <- D2[D2$lysate_reads >0 | D2$homogenate_reads > 0,]

D2_filtered <- D2 %>%
  filter(Phylum %in% "Arthropoda")
length(unique(D2_filtered$cluster))
D2_filtered$taxon <- paste0(D2_filtered$Order, "_", D2_filtered$Family)
D2_filtered$taxon2 <- paste0(D2_filtered$Family, " (", D2_filtered$Order, ")")


#Filter out the NAs for plotting
D2_filtered_plot <- D2_filtered[!is.na(D2_filtered$lys_homog_ratio) & is.finite(D2_filtered$lys_homog_ratio), ]
#Additionally, eliminate all taxa that are unclassifies (which is essentially few Diptera clusters)
D2_filtered_plot <- D2_filtered %>%
  filter(!str_detect(taxon2, "unclassified"))
D2_filtered_plot <- D2_filtered_plot %>%
  filter(!str_detect(taxon2, "_X"))

#double-check that we don't have any lines with 0 reads in both treatments:
D2_filtered_plot_zeros <- D2_filtered_plot %>%
  filter(lysate_reads == 0 | homogenate_reads == 0)

ZERO_summary <- D2_filtered_plot_zeros %>%
  mutate(category = case_when(
    lysate_reads > 0 & homogenate_reads == 0 ~ "lysate_only",
    lysate_reads == 0 & homogenate_reads > 0 ~ "homogenate_only",
    TRUE ~ "zero_in_both"
  )) %>%
  group_by(Family, taxon2, category) %>%
  summarise(n_samples = n(), .groups = "drop")

ZERO_filtered <- ZERO_summary %>%
  filter(category != "zero_in_both")


#### Get the number of families presented in one treatment only that were in the Fig2A (Ratio)
# Families shown in the ratio plot
families_in_p1 <- unique(D_filtered_plot_500$taxon2)

# Zero-only counts, filtered to those families
ZERO_summary <- D2_filtered_plot_zeros %>%
  mutate(category = case_when(
    lysate_reads > 0 & homogenate_reads == 0 ~ "lysate_only",
    lysate_reads == 0 & homogenate_reads > 0 ~ "homogenate_only",
    TRUE ~ "zero_in_both"
  )) %>%
  filter(category != "zero_in_both") %>%
  group_by(Order, Family, taxon2, category) %>%
  summarise(n_clusters = n_distinct(cluster), .groups = "drop") %>%
  filter(taxon2 %in% families_in_p1)   # keep only families in p1

p2 <- ggplot(ZERO_summary, aes(x = n_clusters, y = taxon2, fill = category)) +
  geom_col(position = "dodge") +
  facet_grid(rows = vars(Order), scales = "free_y", space = "free_y") +
  theme_minimal() +
  theme(
    strip.text.y = element_blank(),   # don’t duplicate y labels
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "gray90", color = "black")
  ) +
  scale_fill_manual(values = c("lysate_only" = "#FD6467", "homogenate_only" = "darkgrey")) +
  labs(x = "Number of clusters", y = "", fill = "Category")
p2
ggsave(file="Fig2_treatment_exclusive_clusters.jpg", height=13, width=7, plot = p2)

####
####
####
####


#### for the selected 15 samples only:
SelSamples <- c("S8APX1","SLWUQ2","S5MEZA","SY34ZS","STHFJP","SKYQAL","SIRUBH","SCNXZK","SUXEKB","SQVR9H","SIUQ3G","SG3CB2","SBWPUY","SC1YIM","SACBBR")
D_15Samples <- D %>% filter(sample %in% SelSamples)
unique(D_15Samples$cluster)

# Look at reads per sample for homogs and lysates:
# Sum reads per sample and treatment
read_summary <- D_15Samples %>%
  group_by(sample) %>%
  summarize(
    lysate_total = sum(lysate_reads, na.rm = TRUE),
    homogenate_total = sum(homogenate_reads, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(lysate_total, homogenate_total),
               names_to = "treatment",
               values_to = "total_reads") %>%
  mutate(treatment = ifelse(treatment == "lysate_total", "Lysate", "Homogenate"))

# data inspection - barplot with biomass 
#add biomass data to the plot:
#get biomass data and add it ot the reads summary:
meta<-read.delim2("./15_samples_metadata.csv", sep=";")
read_summary_biomass <- read_summary %>%
  left_join(meta, by = "sample")
read_summary_biomass$Biomass<-as.numeric(read_summary_biomass$Biomass)
# Rescale biomass so it fits visually on read count scale
max_reads <- max(read_summary_biomass$total_reads, na.rm = TRUE)
max_biomass <- max(read_summary_biomass$Biomass, na.rm = TRUE)
read_summary_biomass <- read_summary_biomass %>%
  mutate(biomass_scaled = Biomass * max_reads / max_biomass)
#Plot reads with overlay of biomass:
mass<-ggplot(read_summary_biomass, aes(x = sample)) +
  geom_col(aes(y = total_reads, fill = treatment), position = "dodge") +
  geom_point(aes(y = biomass_scaled), color = "black", shape = 21, size = 2, 
             position = position_dodge(width = 0.9)) +
  scale_y_continuous(
    name = "Total Reads",
    sec.axis = sec_axis(~ . * max_biomass / max_reads, name = "Biomass (mg)")
  ) +
  theme_minimal() +
  labs(
    x = "Sample",
    fill = "Treatment",
    title = "Total Reads per Sample with Biomass"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
mass
ggsave(file="Fig_S5_15samples_reads_biomass.jpg", height=7, width=8, plot = mass)

## Exploring 15 sample data by scatter plots:
## Prep data - filter out anything but Arthropoda
D_15Samples_filtered <- D_15Samples %>%
  filter(Phylum %in% "Arthropoda")
unique(D_15Samples_filtered$cluster)
D_15Samples_filtered$taxon <- paste0(D_15Samples_filtered$Order, "_", D_15Samples_filtered$Family)
D_15Samples_filtered$taxon2 <- paste0(D_15Samples_filtered$Family, " (", D_15Samples_filtered$Order, ")")


plot<-ggplot(D_15Samples_filtered, aes(x = lysate_reads, y = homogenate_reads, color = taxon)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Lysate Reads",
    y = "Homogenate Reads",
    title = "15 samples only",
    color = "Taxonomy"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
plot

legend<-ggplot(D_15Samples_filtered, aes(x = lysate_reads, y = homogenate_reads, color = taxon)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Lysate Reads",
    y = "Homogenate Reads",
    title = "Scatter Plot of Lysate vs. Homogenate Reads",
    color = "Taxonomy"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

#ggsave(file="Fig_matched_15Samples.jpg", height=7, width=8, plot = plot)
#ggsave(file="Fig_matched_15Samples_legend.jpg", height=12, width=15, plot = legend)

f<-ggplot(D_15Samples_filtered, aes(x = lysate_reads, y = homogenate_reads, color = taxon)) +
  geom_point(alpha = 0.7) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ Order) +
  labs(
    x = "Lysate Reads",
    y = "Homogenate Reads",
    title = "Scatter Plot of Lysate vs. Homogenate Reads",
    color = "Taxonomy"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
f

#ggsave(file="Fig_matched_15Samples_legend.jpg", height=12, width=15, plot = f)
#ggsave(file="Fig_matched_15Samples_byORDER.jpg", height=7, width=8, plot = f)


########## Selected 15 samples PLOT a ratio. 
# average read depth varies between lysate and homogenate. find out by how much for 15 samples:
lys15_sample_mean <- sum(D_15Samples$lysate_reads)/length(unique(D_15Samples$sample))
hom15_sample_mean <- sum(D_15Samples$homogenate_reads)/length(unique(D_15Samples$sample))
factor<-lys15_sample_mean/hom15_sample_mean

D_15Samples$homogenate_reads_corrected <- D_15Samples$homogenate_reads * factor
#checking id the means of "lysate reads" and "homogenate_reads_corrected" are the same.
summary(D_15Samples)
#they are. we proceed

## Exclude all pairs of samples that have 0 in any treatment.
#and calculate the ratio of lysate reads to (corrected) homogenate reads.
D_15Samples$lys_homog_ratio <- ifelse(
  D_15Samples$lysate_reads != 0 & D_15Samples$homogenate_reads != 0,
  D_15Samples$lysate_reads / D_15Samples$homogenate_reads_corrected,
  NA  # Use NA for rows where one or both values are zero
)

#Filter anything that is not arthropoda
D_15Samples_filtered <- D_15Samples %>%
  filter(Phylum %in% "Arthropoda")
unique(D_15Samples_filtered$cluster)
D_15Samples_filtered$taxon <- paste0(D_15Samples_filtered$Order, "_", D_15Samples_filtered$Family)
D_15Samples_filtered$taxon2 <- paste0(D_15Samples_filtered$Family, " (", D_15Samples_filtered$Order, ")")


#Filter out the NAs for plotting
D_15Samples_filtered_plot <- D_15Samples_filtered[!is.na(D_15Samples_filtered$lys_homog_ratio) & is.finite(D_15Samples_filtered$lys_homog_ratio), ]

# Filter taxa with at least 10 data points
D_15Samples_filtered_plot_10 <- D_15Samples_filtered_plot %>%
  group_by(taxon) %>%
  filter(n() >= 10) %>%
  ungroup()
#Additionally, eliminate all taxa that are unclassifies (which is essentially few Diptera clusters)
D_15Samples_filtered_plot_10 <- D_15Samples_filtered_plot_10 %>%
  filter(!str_detect(taxon2, "unclassified"))

# Boxplot split by Order
bh<-ggplot(D_15Samples_filtered_plot_10, aes(x = lys_homog_ratio, y = Family, fill = Order)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_log10() +
  facet_grid(rows = vars(Order), scales = "free_y", space = "free_y") +
  theme_minimal() +
  labs(
    x = "Lysate / Homogenate Read Ratio (log scale)",
    y = "",
    title = ""
  ) +
  theme(
    legend.position = "none",
    strip.text.y = element_text(angle = 0)  # Optional: rotate facet strip text if needed
  ) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "gray90", color = "black")
  )
bh
ggsave(file="Fig_S6_15Samples_Top10_Ratio.jpg", height=10, width=7, plot = bh)


## Violin plots:
# define plotting functions
plot_lys <- function(D,lbls) {
    ggplot(data=D, aes(x=cluster, y=log10(lysate_reads+1))) +
        theme_minimal() +
        geom_jitter(color="black",size=0.4,alpha=0.5) +
#        geom_dotplot(binaxis="y", stackdir="center", dotsize=0.4, alpha=0.5, binwidth=0.1) +
        geom_violin(alpha=0.5) +
        ylim(c(-0.05,6.5)) +
        scale_x_discrete(labels=lbls)
}

plot_hom <- function(D,lbls) {
    ggplot(data=D, aes(x=cluster, y=log10(homogenate_reads+1))) +
        theme_minimal() +
        geom_jitter(color="black",size=0.4,alpha=0.5) +
#        geom_dotplot(binaxis="y", stackdir="center", dotsize=0.4, alpha=0.5, binwidth=0.1) +
        geom_violin(alpha=0.5) +
        ylim(c(-0.05,6.5)) +
        scale_x_discrete(labels=lbls)
}

# Keep only the 20 most frequent clusters
x <- sort(table(D$cluster),decreasing=TRUE)
E <- D[D$cluster %in% names(x[1:20]),]
E$cluster <- factor(E$cluster,levels=names(x[1:20]))

# Define x axis labels
x_labels <- substr(names(x[1:20]),1,3)

p1 <- plot_lys(E, x_labels)
p2 <- plot_hom(E, x_labels)
#ggsave(file="Fig_matched_most_occurring.jpg", height=7, width=14, plot = p1 + p2 + plot_layout(ncol=1) + plot_annotation(tag_levels="A"))
p1
p2
# Keep the clusters occurring in more than 100 samples
Fd <- D[D$cluster %in% names(x[x>100]),]
Fd$cluster <- factor(Fd$cluster,levels=names(x[x>100]))

# Aggregate data for clusters so we can find meaningful subgroups

# Mean reads is larger for lysates than homogenates; find out by how much
lys_sample_mean <- sum(D$lysate_reads)/length(unique(D$sample))
hom_sample_mean <- sum(D$homogenate_reads)/length(unique(D$sample))

L <- aggregate(Fd$lysate_reads[Fd$lysate_reads!=0]~Fd$cluster[Fd$lysate_reads!=0],FUN=mean)
H <- aggregate(Fd$homogenate_reads[Fd$homogenate_reads!=0]~Fd$cluster[Fd$homogenate_reads!=0],FUN=mean)
#L <- aggregate(F$lysate_reads~F$cluster,FUN=mean)
#H <- aggregate(F$homogenate_reads~F$cluster,FUN=mean)
colnames(L) <- c("cluster","mean_lys_reads")
colnames(H) <- c("cluster","mean_hom_reads")
LH <- merge(L, H, by="cluster")

set1 <- LH$cluster[order(LH$mean_lys_reads/lys_sample_mean+LH$mean_hom_reads/hom_sample_mean, decreasing=TRUE)]
set2 <- LH$cluster[order(LH$mean_lys_reads/LH$mean_hom_reads,decreasing=TRUE)]
set3 <- LH$cluster[order(LH$mean_lys_reads/LH$mean_hom_reads)]
set4 <- LH$cluster[order(LH$mean_lys_reads/lys_sample_mean+LH$mean_hom_reads/hom_sample_mean)]

cutoff <- 20
set1 <- set1[1:cutoff]
set2 <- set2[1:cutoff]
set3 <- set3[1:cutoff]
set4 <- set4[1:cutoff]

lbl_len <- 5
x_labels <- substr(set1,1,lbl_len)
p3 <- plot_lys(D[D$cluster %in% set1,], x_labels)
p4 <- plot_hom(D[D$cluster %in% set1,], x_labels)
#ggsave(file="Fig_matched_most_reads.jpg", height=7, width=14, plot = p3 + p4 + plot_layout(ncol=1) + plot_annotation(tag_levels="A"))
p3
p4
x_labels <- substr(set2,1,lbl_len)
p5 <- plot_lys(D[D$cluster %in% set2,], x_labels)
p6 <- plot_hom(D[D$cluster %in% set2,], x_labels)
ggsave(file="Fig_S2_matched_largest_lysate_surplus.jpg", height=7, width=14, plot = p5 + p6 + plot_layout(ncol=1) + plot_annotation(tag_levels="A"))
p5
p6
x_labels <- substr(set3,1,lbl_len)
p7 <- plot_lys(D[D$cluster %in% set3,], x_labels)
p8 <- plot_hom(D[D$cluster %in% set3,], x_labels)
ggsave(file="Fig_S3_matched_largest_homogenate_surplus.jpg", height=7, width=14, plot = p7 + p8 + plot_layout(ncol=1) + plot_annotation(tag_levels="A"))
p7
p8
x_labels <- substr(set4,1,lbl_len)
p9 <- plot_lys(D[D$cluster %in% set4,], x_labels)
p10 <- plot_hom(D[D$cluster %in% set4,], x_labels)
#ggsave(file="Fig_matched_least_reads.jpg", height=7, width=14, plot = p9 + p10 + plot_layout(ncol=1) + plot_annotation(tag_levels="A"))
p9
p10
