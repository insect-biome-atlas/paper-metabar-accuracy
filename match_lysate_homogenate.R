# Script for generating table matching lysate and homogenate
# reads for Swedish IBA data
library(data.table)
library(reshape2)

# Read in fxn for getting data
# TODO: Change directory to your local copy of the IBA utils repo
cur_dir <- getwd()
setwd("~/dev/ms-repos-iba/utils/")
source("get_iba_co1_data_fxn.R")
setwd(cur_dir)

# Set IBA data paths
# TODO: Change to your local paths to the IBA data
data_path <- "~/dev/figshare-repos/iba/processed_data/v3/"
metadata_path <- "~/dev/figshare-repos/iba/raw_data/v6/"

# Get IBA lysate and homogenate data
# Note that these are data tables and not standard R data frames
# Spike-ins are removed by default
L <- get_iba_co1_data(data_path, metadata_path, country="SE", dataset="lysate", calibrate=TRUE, data.table=FALSE)
H <- get_iba_co1_data(data_path, metadata_path, country="SE", dataset="homogenate", calibrate=TRUE, data.table=FALSE)
L1 <- data.frame(L)
H1 <- data.frame(H)

# Get basic stats on the lysates:
###Number of reds per sample
x <- colnames(L1)
samples <- grepl("_",x) & grepl("P",x)
Sums<-colSums(L1[samples])
###Filter one failed sample with only 1 read
SumsFilt<- subset(Sums, Sums>1000)

###Mean and Sd
summary(Sums)
str(SumsFilt)
mean(Sums)
mean(SumsFilt)
sd(SumsFilt)


# Get sample metadata
M <- read.delim(paste0(metadata_path,"CO1_sequencing_metadata_SE.tsv"))

# Convert the column sample names to sampleID_FIELD
x <- colnames(H1)
samples <- grepl("_",x) & grepl("P",x)
field_samples <- M$sampleID_FIELD[match(x[samples],M$sampleID_NGI)]
x[samples] <- field_samples
colnames(H1) <- x

x <- colnames(L1)
samples <- grepl("_",x) & grepl("P",x)
field_samples <- M$sampleID_FIELD[match(x[samples],M$sampleID_NGI)]
x[samples] <- field_samples
colnames(L1) <- x

# Extract the lysate data columns matching homogenate data
# Remove any homogenate data columns without matching lysate data (if any)
L2 <- L1[ ,colnames(L1) %in% colnames(H1)]
H2 <- H1[ ,colnames(H1) %in% colnames(L2)]
H2 <- H2[ ,colnames(L2)]    # Order columns in the same way

# Remove clusters absent in these samples
first_sample <- 1 + which(colnames(L2)=="representative")
L2_tot <- rowSums(L2[,first_sample:ncol(L2)])
H2_tot <- rowSums(H2[,first_sample:ncol(H2)])
L3 <- L2[L2_tot > 0,]
H3 <- H2[H2_tot > 0,]

# Write these data frames
write.table(L3, "shared_lysates.tsv", row.names=FALSE, sep="\t")
write.table(H3, "shared_homogenates.tsv", row.names=FALSE, sep="\t")

# Melt datasets
id.vars.idx <- 1:(first_sample-1)
L4 <- reshape2::melt(L3, id.vars=1:(first_sample-1), variable.name="sample", value.name="lysate_reads")
H4 <- reshape2::melt(H3, id.vars=1:(first_sample-1), variable.name="sample", value.name="homogenate_reads")

# Merge datasets
LH <- merge(L4, H4, all=TRUE)
LH$lysate_reads[is.na(LH$lysate_reads)] <- 0
LH$homogenate_reads[is.na(LH$homogenate_reads)] <- 0

# Remove all-zero entries
LH1 <- LH[LH$lysate_reads > 0 | LH$homogenate_reads > 0,]

# Write merged data frame
write.table(LH1, "matched_lysate_homogenate.tsv", row.names=FALSE, sep="\t")


