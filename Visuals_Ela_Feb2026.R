setwd("~/Desktop/git/paper-metabar-accuracy/")

install.packages('ggplot2')
library('ggplot2')
# VEnn of all methods:
install.packages("VennDiagram")
library(VennDiagram)
install.packages("tidyverse")
library("tidyverse")
install.packages("ggvenn")
library("ggvenn") 
install.packages("wesanderson")
library("wesanderson")
install.packages("vegan")
library("vegan")
install.packages("reshape2")
library(reshape2)
library(ggplot2)
library(dplyr)


# Look at barcoding success sample by sample:
Tab_plus_failed=read.delim("Barcoding_cleaned_matched_corrected_plusFAILED.csv", header=TRUE, sep=";", stringsAsFactors=FALSE)

# Create a new category combining Outcome and MatchStatus
df <- Tab_plus_failed %>%
  mutate(Category = case_when(
    Barcoding_outcome == "Failed_barcoding" ~ "Failed_barcoding",
    Barcoding_outcome == "Successful_barcoding" & Matched_to_Cluster == "YES" ~ "Success-Matched",
    Barcoding_outcome == "Successful_barcoding" & Matched_to_Cluster == "NO" ~ "Success-Unmatched"
  ))

# Summarize occurrences per SampleID and Category
summary_df <- df %>%
  count(sample_name, Category)

summary_df2 <- df %>%
  count(Category)

#How many unique barcodes do we have:
length(unique(df$barcode))

df_unmatched <- subset(df, df$Barcoding_outcome == "Success-Unmatched")
length(unique(df_unmatched$barcode))

# Create color palette
colors <- c("#D67236" , "#FD6467","#EFC000FF") # "#5B1A18", 

# Plot histogram
ggplot(summary_df, aes(x = sample_name, y = n, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Barcoding and matching success per sample",
       x = "Sample ID",
       y = "Number of individuals",
       fill = "") +
  scale_fill_manual(values = colors) +  # Set custom colors
  theme_minimal() +  # Minimal background
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x labels for readability

#_____________________________________________________________________________
# Reading in the metabarcoding - clustered-cleaned counts:
# Cluster counts for 45 sampples and 5 controls:
MetaBar_raw <- read.table("cleaned_nochimera_cluster_counts_ELA001_no0_QR.txt", row.names =1, sep = "\t", header=T)
#Metadata file
MetaDat <- read.delim("Metadata_barcoding.txt")
MetaSample <- read.csv("Sample_metadata.csv", sep = ";")
# Filtering of counts file 
# Eliminating positive and negative controls from the counts data file:
MetaBar_samples_raw <- MetaBar_raw[,subset(MetaDat$Sample_ID_QR, MetaDat$Type=="Sample")]
#Getting rid of rows with 0 counts.
MetaBar_samples <-MetaBar_samples_raw[rowSums(MetaBar_samples_raw) > 0, ]

## getting stats on the number of clusters in different treatments.
head(MetaBar_samples)
summary(MetaBar_samples)

### Total number of clusters found in 15 samples metabarcoding (excluding positive and negative cntrls)
nrow(MetaBar_samples)

# Identify treatment columns
lysate_cols <- grepl("_L$", colnames(MetaBar_samples))
homogenate_cols <- grepl("_H$", colnames(MetaBar_samples))

# Function to compute all stats for a treatment
get_treatment_stats <- function(mat) {
  
  # Total reads per sample
  reads_per_sample <- colSums(mat)
  
  # Total clusters present in treatment (at least one non-zero)
  total_clusters <- sum(rowSums(mat) > 0)
  
  c(
    n_samples = ncol(mat),
    mean_reads = mean(reads_per_sample),
    median_reads = median(reads_per_sample),
    sd_reads = sd(reads_per_sample),
    min_reads = min(reads_per_sample),
    max_reads = max(reads_per_sample),
    total_clusters = total_clusters
  )
}

# Apply to both treatments
final_stats <- rbind(
  Lysates = get_treatment_stats(MetaBar_samples[, lysate_cols]),
  Homogenates = get_treatment_stats(MetaBar_samples[, homogenate_cols])
)

final_stats


##############################################
#get lists of clusters that were detected in each of the treatments (lysate, homog)
#### Lysate
lysate_samples <- MetaDat %>%
  filter(Treatment == "Lysis") %>%
  pull(Sample_ID_QR)
# Subset the counts table to only those samples
counts_lysate <- MetaBar_samples[, lysate_samples]
# Find clusters that are detected (non-zero) in at least one of those samples
clusters_detected_lysate <- rownames(counts_lysate)[rowSums(counts_lysate > 0) > 0]
#### Homogenate
homog_samples <- MetaDat %>%
  filter(Treatment == "Homogen") %>%
  pull(Sample_ID_QR)
# Subset the counts table to only those samples
counts_homog <- MetaBar_samples[, homog_samples]
# Find clusters that are detected (non-zero) in at least one of those samples
clusters_detected_homog <- rownames(counts_homog)[rowSums(counts_homog > 0) > 0]

#&&&&&&&&&&&&&&&&&&&&&&&&&
# How many clusters from each treatment (lysis/homog) match with the ground truth / barcoding
# Add to the dataframe:
Tab_plus_failed =read.delim("Barcoding_cleaned_matched_corrected_plusFAILED.csv", header=TRUE, sep=";", stringsAsFactors=FALSE)

Tab_plus_failed <- Tab_plus_failed %>%
  mutate(
    #Ethanol_treatment  = ifelse(CLUSTER_ID_otu %in% clusters_detected_ethanol, "YES", "NO"),
    Lysis_treatment = ifelse(CLUSTER_ID_otu %in% clusters_detected_lysate, "YES", "NO"),
    Homog_treatment    = ifelse(CLUSTER_ID_otu %in% clusters_detected_homog, "YES", "NO")
  )

#### Join this with the MetaData file for sample - to get proper sample Id for the plot
Tab_plus_failed <- Tab_plus_failed %>%
  left_join(MetaSample, by = "sample_name")



# Create a new category combining Outcome and MatchStatus
df <- Tab_plus_failed %>%
  mutate(Category = case_when(
    Barcoding_outcome == "Failed_barcoding" ~ "Failed_barcoding",
    Barcoding_outcome == "Successful_barcoding" & Lysis_treatment == "YES" & Homog_treatment == "NO" ~ "Matched-Lysate_only",
    Barcoding_outcome == "Successful_barcoding" & Lysis_treatment == "YES" & Homog_treatment == "YES" ~ "Matched-Lysate&Homog",
    Barcoding_outcome == "Successful_barcoding" & Homog_treatment == "YES" & Lysis_treatment == "NO" ~ "Matched-Homog_only",
    Barcoding_outcome == "Successful_barcoding" & Lysis_treatment == "NO" & Homog_treatment == "NO" ~ "Unmatched",
    Barcoding_outcome == "Successful_barcoding" & Matched_to_Cluster == "NO" ~ "Unmatched"
  ))

df_clean <- df %>%
  filter(!is.na(Category))

df_clean$taxon2 <- paste0("(", df_clean$Order, ") ", df_clean$Family) 

#how many BARCODES we have overall:
length(unique(df_clean$barcode))

df_clean_matched<-subset(df_clean, df_clean$Matched_to_Cluster == "YES")

# Summarize occurrences per Barcoding outcome
df_unmatched <- df_clean %>%
  group_by(Barcoding_outcome) %>%
  summarise(n = n(), .groups = "drop")


#how many are matched to metabarcoding:
unique_matched_seqs <- df %>%
  filter(Matched_to_Cluster == "YES") %>%
  distinct(barcode)
nrow(unique_matched_seqs)



# So all 4 barcodes found uniquely in EtOH treatment are also all 4 in the set of 947 "matched bardes.
# In essence the number of matched Barcodes is 947-4=943
True_unique_matched_seqs <- setdiff(unique_matched_seqs$barcode, query)

# Make a subset of all df_clean data, that matched
df_clean_True <- df_clean %>%
  filter(barcode %in% True_unique_matched_seqs)
#and how many clusters did we match to? 474.
length(unique(df_clean_True$CLUSTER_ID_otu))


#####
#Get lists of clusters that matched barcoding and that did not match barcoding.





# Summarize occurrences per SampleID and Category
df_summary <- df_clean %>%
  group_by(sample_ID_graph, Category) %>%
  summarise(n = n(), .groups = "drop")


# Create color palette
colors <- c("#D67236" , "darkgrey",  "#FD6467","#EFC000FF", "#5B1A18")

# Plot histogram
ggplot(df_summary, aes(x = sample_ID_graph, y = n, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Barcoding and matching success per sample",
       x = "Sample ID",
       y = "Number of individuals",
       fill = "") +
  scale_fill_manual(values = colors) +  # Set custom colors
  theme_minimal() +  # Minimal background
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x labels for readability


category_summary <- df_clean %>%
  group_by(Category) %>%
  summarise(total = n(), .groups = "drop") %>%
  arrange(desc(total))

category_summary

ggplot(category_summary, aes(x = reorder(Category, -total), y = total, fill = Category)) +
  geom_col() +
  geom_text(aes(label = total), vjust = -0.5) +
  theme_minimal() +
  scale_fill_manual(values = colors) +
  labs(
    title = "Barcoding and matching to metabarcoding results outcomes across all samples",
    x = "",
    y = "Number of individuals"
  ) +
  theme(
    axis.text.x = element_blank(),      # FIXED
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )



# Filter for the three categories of interest
selected_taxa <- df_clean %>%
  filter(Category %in% c("Matched-Lysate_only", "Matched-Homog_only", "Unmatched"))

# Summarise at chosen taxonomic level (e.g., Order–Species)
summary_taxa <- selected_taxa %>%
  group_by(Category, Order, Family, Genus, Species) %>%
  summarise(n_clusters = n_distinct(CLUSTER_ID_otu),
            .groups = "drop")

# Plot: stacked bar per category
ggplot(summary_taxa, aes(x = Category, y = n_clusters, fill = Order)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Taxonomic composition of selected categories",
       x = "Category",
       y = "Number of clusters",
       fill = "Order") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot UNMATCHED clusters only, per sample:
# Filter for 'Unmatched' only
unmatched_taxa <- df_clean %>%
  filter(Category %in% c("Unmatched"))

unique(unmatched_taxa$otu)
unique(unmatched_taxa$barcode)
unique(unmatched_taxa$Order)

write_delim(unmatched_taxa, "List_unmatched_barcodes.txt")

# Summarize 
df_unmatched_summ <- unmatched_taxa %>%
  group_by(sample_ID_graph, Family, Order) %>%
  summarise(n = n(), .groups = "drop")

family_counts <- df_unmatched_summ %>%
  group_by(Order, Family) %>%
  summarise(Total = sum(n), .groups = "drop")

# Sort families within each order
family_counts <- family_counts %>%
  arrange(Order, desc(Total)) %>%
  mutate(Family = factor(Family, levels = unique(Family)))

# Plot: all families, grouped visually by order using facets
un<-ggplot(family_counts, aes(y = Family, x = Total, fill = Order)) +
  geom_col() +
  facet_grid(
    Order ~ ., 
    scales = "free_y", 
    space = "free_y", 
    switch = "y",
    labeller = labeller(Order = label_wrap_gen(width = 20, multi_line = FALSE))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.y.left = element_text(angle = 0, hjust = 0),  # <-- makes order names horizontal
    strip.placement = "outside",
    panel.spacing = unit(0.2, "lines"),
    legend.position = "none"
  ) +
  labs(
    title = "Barcodes without a match to metabarcoding",
    x = "Number of Specimens",
    y = NULL
  )
un
ggsave(filename = "Unmatched_occurences_by_family.jpeg", device="jpeg", width = 7, height = 9, un)




# Plot MATCHED-LYSATE=ONLY clusters only, per sample:
# Filter for 'Unmatched' only
Lysate_taxa <- df_clean %>%
  filter(Category %in% c("Matched-Lysate_only"))

# Summarize 
df_Lysate_summ <- Lysate_taxa %>%
  group_by(sample_ID_graph, Family, Order) %>%
  summarise(n = n(), .groups = "drop")

family_counts_lys <- df_Lysate_summ %>%
  group_by(Order, Family) %>%
  summarise(Total = sum(n), .groups = "drop")

# Sort families within each order
family_counts_lys <- family_counts_lys %>%
  arrange(Order, desc(Total)) %>%
  mutate(Family = factor(Family, levels = unique(Family)))

# Plot: all families, grouped visually by order using facets
ux<-ggplot(family_counts_lys, aes(y = Family, x = Total, fill = Order)) +
  geom_col() +
  facet_grid(
    Order ~ ., 
    scales = "free_y", 
    space = "free_y", 
    switch = "y",
    labeller = labeller(Order = label_wrap_gen(width = 20, multi_line = FALSE))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.y.left = element_text(angle = 0, hjust = 0),  # <-- makes order names horizontal
    strip.placement = "outside",
    panel.spacing = unit(0.2, "lines"),
    legend.position = "none"
  ) +
  labs(
    title = "Barcodes matched to lysate metabarcoding",
    x = "Number of Specimens",
    y = NULL
  )
ux
ggsave(filename = "Lysate_only_occurences_by_family.jpeg", device="jpeg", width = 7, height = 9, ux)


# Plot MATCHED-Homog_ONLY clusters only, per sample:
# Filter for 'Unmatched' only
Hom_taxa <- df_clean %>%
  filter(Category %in% c("Matched-Homog_only"))

# Summarize 
df_Hom_summ <- Hom_taxa %>%
  group_by(sample_ID_graph, Family, Order) %>%
  summarise(n = n(), .groups = "drop")

family_counts_hom <- df_Hom_summ %>%
  group_by(Order, Family) %>%
  summarise(Total = sum(n), .groups = "drop")

# Sort families within each order
family_counts_hom <- family_counts_hom %>%
  arrange(Order, desc(Total)) %>%
  mutate(Family = factor(Family, levels = unique(Family)))

# Plot: all families, grouped visually by order using facets
uh<-ggplot(family_counts_hom, aes(y = Family, x = Total, fill = Order)) +
  geom_col() +
  facet_grid(
    Order ~ ., 
    scales = "free_y", 
    space = "free_y", 
    switch = "y",
    labeller = labeller(Order = label_wrap_gen(width = 20, multi_line = FALSE))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.y.left = element_text(angle = 0, hjust = 0),  # <-- makes order names horizontal
    strip.placement = "outside",
    panel.spacing = unit(0.2, "lines"),
    legend.position = "none"
  ) +
  labs(
    title = "Barcodes matched to homogenate metabarcoding",
    x = "Number of Specimens",
    y = NULL
  )
uh
ggsave(filename = "Homog_only_occurences_by_family.jpeg", device="jpeg", width = 7, height = 9, uh)






# Reading in the barcoding results:
Tab=read.delim("Barcoding_cleaned_matched_corrected.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE, row.names = 2 )
# Number of species and unmatched (i.e. clusterID-identified&matched species- or zOTU-unmatched)
length(unique(Tab$CLUSTER_ID_zotu))
# Number of species and unmatched&clustered (i.e. clusterID-identified&matched species- or OTU-unmatched)
length(unique(Tab$CLUSTER_ID_otu))
# Number of unique sequences
length(unique(Tab$barcode))
# Number of unique zOTUs
length(unique(Tab$zotu))
length(unique(Tab$otu))

#spikes?
#x<-subset(OnlyMatched, Tab$SPIKE_IN == "YES")


# Stats for those barcodes that got a match in metabarcoding. 
OnlyMatched<-subset(Tab, Tab$Matched_to_Cluster == "YES")
length(unique(OnlyMatched$CLUSTER_ID_otu))
length(unique(OnlyMatched$otu))
#spikes?
x<-subset(OnlyMatched, Tab$SPIKE_IN == "YES")


x<-subset(OnlyMatched, OnlyMatched$SPIKE_IN == "YES")
#same number of spiek ins as before - good. cause would have to add up.

#Making a long-format counts for all clusters and OTUs (barcodes that did not get a match in metabarcoding dataset)
counts <- table(Tab$sample_name, Tab$CLUSTER_ID_otu)
c <- data.frame(counts)
c
length(unique(c$Var2))
write.table(c,"Ground_truth_with_UNmatched_OTUs.tsv", row.names = T, sep="\t")

# eliminating all zOTUS without matches.
cd<-read.delim("Ground_truth.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
length(unique(cd$Var2))


# Create a sample data frame
df <- data.frame(Existing_Column = c(3, 8, 2, 6, 4, 10, 1, 9, 3, 7))

unique(Tab$sample_name)
# Define a named vector with values and their corresponding labels
value_labels <- c("SUXEKB" = "Forest_3",
                  "SG3CB2" = "Forest_4",
                  "SQVR9H" = "Forest_1",
                  "SY34ZS" = "Grassland_3",
                  "SBWPUY" = "Forest_2",
                  "SIUQ3G" = "Forest_5",
                  "SKYQAL" = "Grassland_5",
                  "SLWUQ2" = "Wetland_3",
                  "S5MEZA" = "Grassland_2",
                  "SIRUBH" = "Wetland_4",
                  "S8APX1" = "Wetland_2",
                  "SACBBR" = "Wetland_5",
                  "SC1YIM" = "Wetland_1",
                  "STHFJP" = "Grassland_1",
                  "SCNXZK" = "Grassland_4")

# Add a new column based on values in the existing column - labels
Tab$Label <- value_labels[as.character(Tab$sample_name)]


#_____________________________________________________________________________
# Barplots:
#BARCODING only
# Making a graph at order level - Only INSECTA, Spike0ins excluded:
colA2 <- c("#D67236", "#ECCBAE", "#FD6467", "#5B1A18", "#F8AFA8", "#9A8822", "#046C9A", "#ABDDDE", "#74A089", "#D69C4E", "#D8A499", "#7294D4", "#000000", "#C6CDF7", "#DD8D29",  "#46ACC8","#E2D200",  "#E58601", "#B40F20", "#E6A0C4", "#F5CDB4", "#FDDDA0")

matched <- subset(Tab, Tab$Matched_to_Cluster == "YES" & Tab$Class=="Insecta" & Tab$SPIKE_IN=="NO")
CouC <- table(matched$sample_name, matched$Order)
d <- data.frame(CouC)
ggplot(d, aes(x = Var1, y=Freq, fill = Var2)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual("Order", values = colA2) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Who's in there?")

# Other classes (non-insecta) included
colA3 <- c("#000000", "#D67236", "#ECCBAE", "#FD6467", "#B40F20", "#5B1A18", "#F8AFA8", "#9A8822", "#E2D200", "#046C9A", "#FDDDA0", "#ABDDDE", "#7294D4","#74A089",  "#D69C4E",  "#D8A499","#C6CDF7", "#DD8D29",  "#46ACC8", "#E58601",  "#E6A0C4", "#F5CDB4")
matched2 <- subset(Tab, Tab$Matched_to_Cluster == "YES" & Tab$SPIKE_IN=="NO")
CouC2 <- table(matched2$sample_name, matched2$Order)
d2 <- data.frame(CouC2)
ggplot(d2, aes(x = Var1, y=Freq, fill = Var2)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual("Order", values = colA3) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Who's in there?")


colA3 <- c("#C6CDF7", "#000000", "#D67236", "#ECCBAE", "#46ACC8", "#FD6467", "#B40F20", "#5B1A18", "#F8AFA8", "#9A8822", "#E2D200", "#046C9A", "#FDDDA0", "#ABDDDE","#E58601", "#7294D4","#E6A0C4", "#74A089", "#DD8D29", "#F5CDB4", "#D69C4E",  "#D8A499",  "#46ACC8",   "#E6A0C4", "indianred3")
NoSpike <- subset(Tab, Tab$SPIKE_IN=="NO")
CouC3 <- table(NoSpike$Label, NoSpike$Order)
d3 <- data.frame(CouC3)
d3p<-ggplot(d3, aes(x = Var1, y=Freq, fill = Var2)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual("Order", values = colA3) +
  labs(title = "Who's in there? (orders)", x = NULL, y = "Number of individuals") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "Whos_in_there_all_orders.pdf", device="pdf", width = 8, height = 5,  d3p)
ggsave(filename = "Whos_in_there_all_orders.jpeg", device="jpeg", width = 8, height = 5,  d3p)



#_____________________________________________________________________________
# Reading in the metabarcoding - clustered-cleaned counts:
# Cluster counts for 45 sampples and 5 controls:
MetaBar_raw <- read.table("cleaned_nochimera_cluster_counts_ELA001_no0_QR.txt", row.names =1, sep = "\t", header=T)

#Metadata file
MetaDat <- read.delim("Metadata_barcoding.txt")

# Not sure I need this paragraph now.
# Cluster counts for only the clusters that were matched to barcodes:
#  MetaBar_MATCHEDonly <- read.table("cleaned_nochimera_MATCHED_cluster_counts_QRnames.txt", row.names =1, sep = "\t", header=T)

# Filtering of counts file 
# Eliminating positive and negative controls from the counts data file:
MetaBar_samples_raw <- MetaBar_raw[,subset(MetaDat$Sample_ID_QR, MetaDat$Type=="Sample")]
#Getting rid of rows with 0 counts.
MetaBar_samples <-MetaBar_samples_raw[rowSums(MetaBar_samples_raw) > 0, ]

## Filter the Ethanol samples as well! remove all columns ending with "_E"
MetaBar_sample <- MetaBar_samples[, !grepl("_E$", colnames(MetaBar_samples))]


##########----------------------
# VENN DIAGRAMS

# Venn 1) Barcodimg vs. Metabarcodig (lumping all treatments and samples together):

# Get all clusters names:
MetaBarClusts<-row.names(MetaBar_sample)
#Get all matched and unmatched cluster names form Barcoding
BarcClusts<-unique(Tab$CLUSTER_ID_otu)

# Create the Venn diagram
## GGVENN might be my preferred method actually:
W <- list("Barcoding"= BarcClusts, "Metabarcoding"=MetaBarClusts)

v<-ggvenn(W, fill_color = c("#FD6467", "#9A8822"),
          stroke_size = 0.7, set_name_size = 4) +
  ggtitle(paste("All samples"))
ggsave(filename = paste0("Venn_Barcoding_vs_Metabarcoding_all_samples_SMALLER.pdf"), device="pdf", width = 5, height = 3,  v)
v
# move to new plotting page 
grid.newpage() 

# Venn 2) Barcoding vs Three different Metabarcoding methods:
# subset counts table, filter dead-rows (with only 0s) and get the rownames for each condition.
# ethanol
# Lysates
Lysate_samples <- subset(MetaDat$Sample_ID_QR, MetaDat$Treatment == "Lysis")
MetaBar_All_Lysis <- rownames(MetaBar_sample)[rowSums(MetaBar_sample[,Lysate_samples]) != 0]
# homogenates
Homog_samples <- subset(MetaDat$Sample_ID_QR, MetaDat$Treatment == "Homogen")
MetaBar_All_Homog <- rownames(MetaBar_sample)[rowSums(MetaBar_sample[,Homog_samples]) != 0]

BarcClusts           <- unique(BarcClusts)
MetaBar_All_Lysis    <- unique(MetaBar_All_Lysis)
MetaBar_All_Homog    <- unique(MetaBar_All_Homog)
BarcClusts

# drawing Venn
C <-  list("Individual Barcoding"=BarcClusts, "Metabar. Lysates"=MetaBar_All_Lysis, "Metabar. Homogenates"=MetaBar_All_Homog)
cf<- ggvenn(C, fill_color = c("#5B1A18", "#FD6467", "#EFC000FF"),   #"#EFC000FF"
            stroke_size = 0.5, set_name_size = 4) +
  ggtitle(paste("Meta & Barcoding - All 15 samples"))
ggsave(filename = paste0("Venn_methods_all_samples.pdf"), device="pdf", width = 8, height = 5,  cf)
cf
ggsave(filename = paste0("Venn_methods_all_samples.jpg"), device="jpg", width = 8, height = 5,  cf)

### Get lists of clusters specific to only one method:
Barcoding_only <- setdiff(
  BarcClusts,
  union(MetaBar_All_Lysis, MetaBar_All_Homog)
)

Lysis_only <- setdiff(
  MetaBar_All_Lysis,
  union(BarcClusts, MetaBar_All_Homog)
)

Homog_only <- setdiff(
  MetaBar_All_Homog,
  union(BarcClusts, MetaBar_All_Lysis)
)

#and the set overlap between lysis and homogenization:
Lysis_Homog <- setdiff(
  intersect(MetaBar_All_Lysis, MetaBar_All_Homog),
  BarcClusts
)
# overlap between all.
B_L_H   <- Reduce(intersect, list(BarcClusts, MetaBar_All_Lysis, MetaBar_All_Homog))
#overlap Lysis and Homoges (but not barcoding)
L_H     <- setdiff(intersect(MetaBar_All_Lysis, MetaBar_All_Homog), BarcClusts)

#overlap between lysis and barcoding but not in homogenization:
Lysis_Barc <- setdiff(
  intersect(MetaBar_All_Lysis, BarcClusts),
  MetaBar_All_Homog
)

#overlap between homogenization and barcoding but not in lysis :
Hom_Barc <- setdiff(
  intersect(MetaBar_All_Homog, BarcClusts),
  MetaBar_All_Lysis
)

# Prepare for plot as stacked barplot:
region_counts <- data.frame(
  Subset = c("Only Lysis",
             "Only Homogenate",
             "Lysis + Homogenate",
             "All three",
             "Barcoding + Lysis",
             "Barcoding + Homog.",
             "Only Barcoding"
              ),
  Count = c(length(Lysis_only),
            length(Homog_only),
            length(L_H),
            length(B_L_H),
            length(Lysis_Barc),
            length(Hom_Barc),
            length(Barcoding_only)
             ))

# 1) Force the order you want
region_counts$Subset <- factor(
  region_counts$Subset,
  levels = c("Only Lysis",
             "Only Homogenate",
             "Lysis + Homogenate",
             "All three",
             "Barcoding + Lysis",
             "Barcoding + Homog.",
             "Only Barcoding")
)


bs<-ggplot(region_counts, aes(x = "All methods", y = Count, fill = Subset)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5),
            size = 4) +
  scale_fill_manual(values= c("Only Lysis" = "#FD6467",
                              "Only Homogenate"="darkgrey",
                              "Lysis + Homogenate"="#74A089",
                              "All three"="#EFC000FF",
                              "Barcoding + Lysis" =  "#FD6467",
                              "Barcoding + Homog."= "darkgrey",
                              "Only Barcoding"="#7294D4")) +
  coord_flip() +
  theme_classic() +  # removes background grid
  labs(x = "", y = "") +
  theme_void(
  )
bs
ggsave(filename = "Clusters_subsets_barplot.jpg", device="jpg", width = 8, height = 1,  bs)

  

  ### LOAD IN TAXONOMY FILE FOR ALL METABARCODING CLLUSTERS:
MetaBar_TAX <- read.table("cleaned_nochimera_cluster_taxonomy_ELA001.tsv", row.names =1, sep = "\t", header=F)
head(MetaBar_TAX)

Tax_Lysis_only <- MetaBar_TAX[Lysis_only, 1:5]
colnames(Tax_Lysis_only) <- c("Kingdom", "Phylum", "Class", "Order", "Family")
Tax_Homog_only <- MetaBar_TAX[Homog_only, 1:5]
colnames(Tax_Homog_only) <- c("Kingdom", "Phylum", "Class", "Order", "Family")
Tax_Lysis_Homog <- MetaBar_TAX[Lysis_Homog, 1:5]
colnames(Tax_Lysis_Homog) <- c("Kingdom", "Phylum", "Class", "Order", "Family")


summarise_tax <- function(tax_df, group_name) {
  tax_df %>%
    count(Phylum, Class, Order, Family) %>%   # count number of clusters per taxon
    mutate(Group = group_name)
}

df_Lysis  <- summarise_tax(Tax_Lysis_only, "Lysis_only")
df_Homog  <- summarise_tax(Tax_Homog_only, "Homog_only")
df_Shared <- summarise_tax(Tax_Lysis_Homog, "Lysis_Homog")
df_Shared
Tax_summary <- bind_rows(df_Lysis, df_Homog, df_Shared)


Tax_summary
# Aggregate number of clusters per Order + Family + Group
Tax_plot_df <- Tax_summary %>%
  group_by(Group, Order, Family) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(Order, Family)


# 1 Create merged label
Tax_plot_df <- Tax_plot_df %>%
  mutate(Order_Family = paste(Order, Family, sep = " - "))

# 2 Plot horizontal barplot
ggplot(Tax_plot_df, aes(x = Order_Family, y = n, fill = Order)) +
  geom_col(position = "dodge") +
  facet_wrap(~Group, scales = "free_y", ncol = 1) +
  coord_flip() +  # flip axes: labels on y-axis
  theme_bw() +
  theme(axis.text.y = element_text(hjust = 1)) +
  labs(x = "Number of clusters", y = "Order-Family", fill = "Order")



plot_tax_group <- function(data, group_name) {
  
  df <- data %>%
    filter(Group == group_name) %>%
    arrange(desc(n)) %>%
    mutate(Order_Family = factor(Order_Family, levels = rev(Order_Family)))
  
  ggplot(df, aes(x = Order_Family, y = n, fill = Order)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(
      title = group_name,
      x = "Order-Family",
      y = "Number of clusters",
      fill = "Order"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold")
    )
}

p_lysis  <- plot_tax_group(Tax_plot_df, "Lysis_only")
p_homog  <- plot_tax_group(Tax_plot_df, "Homog_only")
p_shared <- plot_tax_group(Tax_plot_df, "Lysis_Homog")

ggsave(filename = paste0("Histogram_metabar_unmatched_Lysis_only.jpg"), device="jpg", width = 5, height = 6,  p_lysis)
ggsave(filename = paste0("Histogram_metabar_unmatched_Homog_only.jpg"), device="jpg", width = 5, height = 6,  p_homog)
ggsave(filename = paste0("Histogram_metabar_unmatched_Lysis_Homog.jpg"), device="jpg", width = 5, height = 6,  p_shared)


###  Historgram to show all the metaarcoding reads per sample and how much of them are coming from clusters that are matched to ground truth. and how much is not.
# metabarcoding counts table
MetaBar_samples
# Start with Lysates:
MetaBar_samples_L <- MetaBar_samples[ ,
                                      grep("_L$", colnames(MetaBar_samples))
]

# 63 clusters that are lysis only :
Lysis_only
# 109 clusters homogen and lysis
Lysis_Homog
#sanity check
intersect(Lysis_only, rownames(MetaBar_samples_L))
intersect(Lysis_Homog, rownames(MetaBar_samples_L))
all_clusters <- rownames(MetaBar_samples_L)
Matched_clusters <- setdiff(all_clusters,
                            union(Lysis_only, Lysis_Homog))
# get all matched and two sub-groups
reads_Lysis <- colSums(MetaBar_samples_L[Lysis_only, , drop = FALSE])
reads_Shared <- colSums(MetaBar_samples_L[Lysis_Homog, , drop = FALSE])
reads_Matched <- colSums(MetaBar_samples_L[Matched_clusters, , drop = FALSE])

df_reads <- data.frame(
  Sample = names(reads_Lysis),
  Lysis_uniq = reads_Lysis,
  Lysis_Homog = reads_Shared,
  Matched = reads_Matched
)

df_long <- df_reads %>%
  pivot_longer(-Sample,
               names_to = "Cluster_Group",
               values_to = "Reads")

lys<- ggplot(df_long, aes(x = Sample, y = Reads, fill = Cluster_Group)) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(values = c(
    "Lysis_uniq" = "#FD6467",   # warm orange
    "Lysis_Homog" = "#5B1A18",   # blue
    "Matched"  = "#EFC000FF"     # neutral background
  )) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Sample",
       y = "Number of reads",
       fill = "Cluster group")
lys
ggsave(filename = paste0("Histogram_metabar_unmatched_counts_.jpg"), device="jpg", width = 5, height = 6,  lys)

### calulate percentage that "spurious" are compared to the total number of counts:

# 1) Total reads across all samples and groups
total_reads_all <- sum(df_long$Reads)

# 2) Sum reads per cluster group
group_totals <- df_long %>%
  group_by(Cluster_Group) %>%
  summarise(Total_Reads = sum(Reads), .groups = "drop")

# 3) Calculate percentage of total
group_totals <- group_totals %>%
  mutate(Percentage_of_total = 100 * Total_Reads / total_reads_all)

group_totals

# Homogenates:
MetaBar_samples_H <- MetaBar_samples[ ,
                                      grep("_H$", colnames(MetaBar_samples))
]

# 89 clusters that are homog only :
Homog_only
# 109 clusters homogen and lysis
Lysis_Homog
#sanity check
intersect(Homog_only, rownames(MetaBar_samples_H))
intersect(Lysis_Homog, rownames(MetaBar_samples_H))
all_clusters <- rownames(MetaBar_samples_H)
Matched_clusters <- setdiff(all_clusters,
                            union(Homog_only, Lysis_Homog))
# get all matched and two sub-groups
reads_Homog <- colSums(MetaBar_samples_H[Homog_only, , drop = FALSE])
reads_Shared <- colSums(MetaBar_samples_H[Lysis_Homog, , drop = FALSE])
reads_Matched <- colSums(MetaBar_samples_H[Matched_clusters, , drop = FALSE])

df_reads <- data.frame(
  Sample = names(reads_Homog),
  Homog_uniq = reads_Homog,
  Lysis_Homog = reads_Shared,
  Matched = reads_Matched
)

df_long <- df_reads %>%
  pivot_longer(-Sample,
               names_to = "Cluster_Group",
               values_to = "Reads")

hom<-ggplot(df_long, aes(x = Sample, y = Reads, fill = Cluster_Group)) +
  geom_col() +
  theme_bw() +
  scale_fill_manual(values = c(
    "Homog_uniq" = "#FD6467",   # warm orange
    "Lysis_Homog" = "#5B1A18",   # blue
    "Matched"  = "#EFC000FF"     # neutral background
  )) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Sample",
       y = "Number of reads",
       fill = "Cluster group")
hom
ggsave(filename = paste0("Histogram_metabar_unmatched_counts_H.jpg"), device="jpg", width = 5, height = 6,  hom)
###
# 1) Total reads across all samples and groups
total_reads_all <- sum(df_long$Reads)

# 2) Sum reads per cluster group
group_totals <- df_long %>%
  group_by(Cluster_Group) %>%
  summarise(Total_Reads = sum(Reads), .groups = "drop")

# 3) Calculate percentage of total
group_totals <- group_totals %>%
  mutate(Percentage_of_total = 100 * Total_Reads / total_reads_all)

group_totals

#### Plot how many 
# 1) Keep only clusters that exist in the count table
lysis_only_present  <- intersect(Lysis_only, rownames(MetaBar_samples_L))
lysis_homog_present <- intersect(Lysis_Homog, rownames(MetaBar_samples_L))

# 2) Subset count table
counts_lysis_only  <- MetaBar_samples_L[lysis_only_present, , drop = FALSE]
counts_lysis_homog <- MetaBar_samples_L[lysis_homog_present, , drop = FALSE]

# 3) Count number of clusters (>0 reads) per sample
richness_lysis_only  <- colSums(counts_lysis_only > 0)
richness_lysis_homog <- colSums(counts_lysis_homog > 0)

# 4) Create dataframe
df_richness_L <- data.frame(
  Sample = names(richness_lysis_only),
  Lysis_only  = richness_lysis_only,
  Lysis_Homog = richness_lysis_homog
)

# 5) Convert to long format
df_long <- df_richness %>%
  pivot_longer(cols = c(Lysis_only, Lysis_Homog),
               names_to = "Cluster_Set",
               values_to = "N_clusters")

# 6) Optional: sort samples by total richness
df_long <- df_long %>%
  group_by(Sample) %>%
  mutate(total = sum(N_clusters)) %>%
  ungroup() %>%
  arrange(desc(total)) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# 7) Stacked barplot
ggplot(df_long, aes(x = Sample, y = N_clusters, fill = Cluster_Set)) +
  geom_col() +   # default = stacked
  coord_flip() +
  scale_fill_manual(values = c(
    "Lysis_only"  = "#FD6467",
    "Lysis_Homog" = "#5B1A18"
  )) +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Number of clusters",
    fill = "Cluster set"
  )



# 1) Keep only clusters that exist in the count table
homog_only_present  <- intersect(Homog_only, rownames(MetaBar_samples_H))
lysis_homog_present <- intersect(Lysis_Homog, rownames(MetaBar_samples_H))

# 2) Subset count table
counts_homog_only  <- MetaBar_samples_H[homog_only_present, , drop = FALSE]
counts_lysis_homog <- MetaBar_samples_H[lysis_homog_present, , drop = FALSE]

# 3) Count number of clusters (>0 reads) per sample
richness_homog_only  <- colSums(counts_homog_only > 0)
richness_lysis_homog <- colSums(counts_lysis_homog > 0)

# 4) Create dataframe
df_richness_H <- data.frame(
  Sample = names(richness_homog_only),
  Homog_only  = richness_homog_only,
  Lysis_Homog = richness_lysis_homog
)

# 5) Convert to long format
df_long <- df_richness %>%
  pivot_longer(cols = c(Homog_only, Lysis_Homog),
               names_to = "Cluster_Set",
               values_to = "N_clusters")

# 6) Optional: sort samples by total richness
df_long <- df_long %>%
  group_by(Sample) %>%
  mutate(total = sum(N_clusters)) %>%
  ungroup() %>%
  arrange(desc(total)) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# 7) Stacked barplot
ggplot(df_long, aes(x = Sample, y = N_clusters, fill = Cluster_Set)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(
    "Homog_only"  = "#EFC000FF",
    "Lysis_Homog" = "#5B1A18"
  )) +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Number of clusters",
    fill = "Cluster set"
  )



### Getting X-Y plot for number of spiurious clusters reads and numbers of specimens failed in barcoding.
failed_counts <- df_summary %>%
  filter(Category == "Failed_barcoding") %>%
  group_by(sample_ID_graph) %>%
  summarise(Failed_barcoding_n = sum(n), .groups = "drop")

failed_counts 

richness_homog_only
richness_lysis_homog
richness_lysis_only

# Convert named vectors to dataframes
df_homog_only <- data.frame(
  Sample = names(richness_homog_only),
  Homog_only = as.numeric(richness_homog_only)
)

df_lysis_homog <- data.frame(
  Sample = names(richness_lysis_homog),
  Lysis_Homog = as.numeric(richness_lysis_homog)
)

df_lysis_only <- data.frame(
  Sample = names(richness_lysis_only),
  Lysis_only = as.numeric(richness_lysis_only)
)

# Remove suffix (_H or _L) to create base sample ID
df_homog_only  <- df_homog_only  %>% mutate(sample_ID = str_remove(Sample, "_H$"))
df_lysis_homog <- df_lysis_homog %>% mutate(sample_ID = str_remove(Sample, "_H$"))
df_lysis_only  <- df_lysis_only  %>% mutate(sample_ID = str_remove(Sample, "_L$"))

# Merge all together by sample_ID
combined_richness <- df_homog_only %>%
  select(sample_ID, Homog_only) %>%
  full_join(df_lysis_homog %>% select(sample_ID, Lysis_Homog),
            by = "sample_ID") %>%
  full_join(df_lysis_only %>% select(sample_ID, Lysis_only),
            by = "sample_ID") %>%
  mutate(across(everything(), ~replace_na(., 0)))

combined_richness


combined_richness <- combined_richness %>%
  left_join(
    MetaSample %>% select(sample_name, sample_ID_graph),
    by = c("sample_ID" = "sample_name")
  )

combined_richness <- combined_richness %>%
  mutate(
    Total_lysis = Lysis_only + Lysis_Homog,
    Total_homog = Homog_only + Lysis_Homog
  )


combined_richness <- combined_richness %>%
  left_join(failed_counts, by = "sample_ID_graph") %>%
  mutate(Failed_barcoding_n = replace_na(Failed_barcoding_n, 0))

library(tidyr)

df_scatter <- combined_richness %>%
  select(sample_ID_graph,
         Failed_barcoding_n,
         Total_lysis,
         Total_homog) %>%
  pivot_longer(
    cols = c(Total_lysis, Total_homog),
    names_to = "Extraction",
    values_to = "Total_clusters"
  )

r2_df <- df_scatter %>%
  group_by(Extraction) %>%
  summarise(
    R2 = summary(lm(Total_clusters ~ Failed_barcoding_n))$r.squared
  )

r2_df <- r2_df %>%
  mutate(
    label = paste0("  R² = ", round(R2, 3)),
    x = max(df_scatter$Failed_barcoding_n),
    y = c(
      max(df_scatter$Total_clusters),
      max(df_scatter$Total_clusters) * 0.9
    )
  )
sc <- ggplot(df_scatter,
               aes(x = Failed_barcoding_n,
                   y = Total_clusters,
                   color = Extraction)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  scale_color_manual(values = c(
    "Total_lysis" = "#FD6467",
    "Total_homog" = "#5B1A18"
  )) +
  labs(x = "Number of failed barcoding specimens",
       y = "Number of unmatched OTUs (>50 reads)",
       color = "Treatment") +
  geom_text(data = r2_df,
            aes(x = x, y = y, label = label, color = Extraction),
            hjust = 1,
            vjust = 12,
            show.legend = FALSE,
            size = 5)

sc
ggsave(filename = paste0("Scatter_failed_brc_unmatched_OTUs.jpg"), device="jpg", width = 5, height = 3,  sc)



### Getting X-Y plot for number of spiurious clusters BUT WITH COUNT >50 READS PER OCCURENCE and numbers of specimens failed in barcoding.
# Lysis_only clusters
lysis_only_present <- intersect(Lysis_only, rownames(MetaBar_samples_L))
counts_lysis_only <- MetaBar_samples_L[lysis_only_present, , drop = FALSE]

richness_lysis_only_50 <- colSums(counts_lysis_only > 50)

# Lysis_Homog clusters in Lysis samples
lysis_homog_present <- intersect(Lysis_Homog, rownames(MetaBar_samples_L))
counts_lysis_homog_L <- MetaBar_samples_L[lysis_homog_present, , drop = FALSE]

richness_lysis_homog_L_50 <- colSums(counts_lysis_homog_L > 50)

# Homog_only clusters
homog_only_present <- intersect(Homog_only, rownames(MetaBar_samples_H))
counts_homog_only <- MetaBar_samples_H[homog_only_present, , drop = FALSE]

richness_homog_only_50 <- colSums(counts_homog_only > 50)

# Lysis_Homog clusters in Homog samples
counts_lysis_homog_H <- MetaBar_samples_H[lysis_homog_present, , drop = FALSE]

richness_lysis_homog_H_50 <- colSums(counts_lysis_homog_H > 50)

Total_lysis_50  <- richness_lysis_only_50 + richness_lysis_homog_L_50
Total_homog_50  <- richness_homog_only_50 + richness_lysis_homog_H_50

library(dplyr)
library(stringr)

df_lysis_50 <- data.frame(
  sample_ID = str_remove(names(Total_lysis_50), "_L$"),
  Total_lysis_50 = as.numeric(Total_lysis_50)
)

df_homog_50 <- data.frame(
  sample_ID = str_remove(names(Total_homog_50), "_H$"),
  Total_homog_50 = as.numeric(Total_homog_50)
)

combined_richness <- combined_richness %>%
  left_join(df_lysis_50, by = "sample_ID") %>%
  left_join(df_homog_50, by = "sample_ID") %>%
  mutate(across(c(Total_lysis_50, Total_homog_50),
                ~replace_na(., 0)))


library(tidyr)
library(ggplot2)

df_scatter_50 <- combined_richness %>%
  select(sample_ID_graph,
         Failed_barcoding_n,
         Total_lysis_50,
         Total_homog_50) %>%
  pivot_longer(
    cols = c(Total_lysis_50, Total_homog_50),
    names_to = "Extraction",
    values_to = "Total_clusters"
  )

#Get r2
r2_df_50 <- df_scatter_50 %>%
  group_by(Extraction) %>%
  summarise(
    R2 = summary(lm(Total_clusters ~ Failed_barcoding_n))$r.squared
  )

r2_df_50 <- r2_df_50 %>%
  mutate(
    label = paste0("  R² = ", round(R2, 3)),
    x = max(df_scatter_50$Failed_barcoding_n),
    y = c(
      max(df_scatter_50$Total_clusters),
      max(df_scatter_50$Total_clusters) * 0.9
    )
  )
sc50 <- ggplot(df_scatter_50,
               aes(x = Failed_barcoding_n,
                   y = Total_clusters,
                   color = Extraction)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  scale_color_manual(values = c(
    "Total_lysis_50" = "#FD6467",
    "Total_homog_50" = "#5B1A18"
  )) +
  labs(x = "Number of failed barcoding specimens",
       y = "Number of unmatched OTUs (>50 reads)",
       color = "Treatment") +
  geom_text(data = r2_df_50,
            aes(x = x, y = y, label = label, color = Extraction),
            hjust = 1,
            vjust = 2,
            show.legend = FALSE,
            size = 5)

sc50
ggsave(filename = paste0("Scatter_failed_brc_unmatched_OTUs_>50.jpg"), device="jpg", width = 5, height = 3,  sc50)



###
install.packages("ggVennDiagram")
yes
library("ggVennDiagram")

lysis_clusters  <- rownames(MetaBar_samples_L)[rowSums(MetaBar_samples_L > 0) > 0]
homog_clusters  <- rownames(MetaBar_samples_H)[rowSums(MetaBar_samples_H > 0) > 0]

venn_list <- list(
  Lysis = lysis_clusters,
  Homogenate = homog_clusters
)

ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#5B1A18") +
  theme_void()

install.packages("eulerr")
yes
library("eulerr")

fit <- euler(venn_list)

plot(fit,
     fills = c("#FD6467", "#5B1A18"),
     edges = TRUE,
     labels = TRUE)


# Do taxonomy plts for low abundance clusters:
#Filter metabarcoding 
low_lysis_clusters <- rownames(MetaBar_samples_L)[
  apply(MetaBar_samples_L, 1, max) < 50
]

low_homog_clusters <- rownames(MetaBar_samples_H)[
  apply(MetaBar_samples_H, 1, max) < 50
]

low_Lysis_only  <- setdiff(low_lysis_clusters, low_homog_clusters)
low_Homog_only  <- setdiff(low_homog_clusters, low_lysis_clusters)
low_Lysis_Homog <- intersect(low_lysis_clusters, low_homog_clusters)

MetaBar_TAX
MetaBar_TAX$Cluster<-row.names(MetaBar_TAX)
colnames(MetaBar_TAX) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "BOLDbin", "Cluster")

Tax_plot_df
low_tax_df <- bind_rows(
  data.frame(Cluster = low_Lysis_only,  Group = "Lysis_only"),
  data.frame(Cluster = low_Homog_only,  Group = "Homog_only"),
  data.frame(Cluster = low_Lysis_Homog, Group = "Lysis_Homog")
) %>%
  left_join(MetaBar_TAX, by = "Cluster") %>%
  count(Group, Order, Family, name = "n") %>%
  mutate(Order_Family = paste(Order, Family, sep = "_"))

p_lysis_low  <- plot_tax_group(low_tax_df, "Lysis_only")
p_homog_low  <- plot_tax_group(low_tax_df, "Homog_only")
p_shared_low <- plot_tax_group(low_tax_df, "Lysis_Homog")
ggsave(filename = paste0("Histogram_metabar_unmatched_Lysis_only_<50.jpg"), device="jpg", width = 5, height = 6,  p_lysis_low)
ggsave(filename = paste0("Histogram_metabar_unmatched_Homog_only_<50.jpg"), device="jpg", width = 5, height = 6,  p_homog_low)
ggsave(filename = paste0("Histogram_metabar_unmatched_Lysis_Homog_<50.jpg"), device="jpg", width = 7, height = 9,  p_shared_low)

p_lys
p_lysis
p_lysis_low

p_homog_low
p_homog
