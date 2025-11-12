# Generate statistics for matches of lysate and homogenate reads

# Name file
file <- "match_stats.txt"

# Get matched data
D <- read.delim("matched_lysate_homogenate.tsv")

# Just in case the matched data include all-zero entries (there is a lot of those)
D <- D[D$lysate_reads > 0 | D$homogenate_reads > 0,]

# Extract taxonomic info about the clusters
Ts <- D[!duplicated(D$cluster),!grepl("reads",colnames(D))]

# Define function for printing taxonomic summary for a cluster set
tax_stats <- function(spec,clusters) {
    cat("\nOrder stats for",spec,"\n",file=file,append=TRUE)
    x <- sort(table(Ts$Order[Ts$cluster %in% clusters]),decreasing=TRUE)
    if (length(x)>20)
        x <- x[1:20]
    print_table(x)

    cat("\nFamily stats for",spec,"\n",file=file,append=TRUE)
    x <- sort(table(Ts$Family[Ts$cluster %in% clusters]),decreasing=TRUE)
    if (length(x)>20)
        x <- x[1:20]
    print_table(x)

    cat("\nGenus stats for",spec,"\n",file=file,append=TRUE)
    x <- sort(table(Ts$Genus[Ts$cluster %in% clusters]),decreasing=TRUE)
    if (length(x)>20)
        x <- x[1:20]
    print_table(x)
}

print_table <- function(x) {
    for (i in 1:length(x))
        cat(names(x)[i],"--",x[i],"\n",file=file,append=TRUE)
}

# Print number of unique clusters
lysates <- unique(D$cluster[D$lysate_reads > 0])
homogenates <- unique(D$cluster[D$homogenate_reads > 0])

cat("There are",length(lysates),"clusters detected in lysates\n",file=file)
cat("There are",length(homogenates),"clusters detected in homogenates\n",file=file,append=TRUE)

cat("Of these clusters,",sum(homogenates %in% lysates),"clusters are shared\n",file=file,append=TRUE)

cat("There are",sum(!(homogenates %in% lysates)),"homogenate clusters missed in lysates\n",file=file,append=TRUE)
cat("There are",sum(!(lysates %in% homogenates)),"lysate clusters missed in homogenates\n",file=file,append=TRUE)

cat("There are",sum(D$lysate_reads==0),"cluster occurrences detected in homogenates but not lysates\n",file=file,append=TRUE)
cat("This represents",round(100*sum(D$lysate_reads==0)/nrow(D),digits=2),"% of all cluster occurrences\n",file=file,append=TRUE)
cat("There are",sum(D$homogenate_reads==0),"cluster occurrences detected in lysates but not homogenates\n",file=file,append=TRUE)
cat("This represents",round(100*sum(D$homogenate_reads==0)/nrow(D),digits=2),"% of all cluster occurrences\n",file=file,append=TRUE)

uniq_lys <- lysates[!(lysates %in% homogenates)]
uniq_hom <- homogenates[!(homogenates %in% lysates)]

tax_stats("clusters found only in lysates", uniq_lys)
tax_stats("clusters found only in homogenates", uniq_hom)

lys_miss <- unique(D$cluster[D$lysate_reads==0])
hom_miss <- unique(D$cluster[D$homogenate_reads==0])

tax_stats("cluster occurrences missed in lysates",lys_miss)
tax_stats("cluster occurrences missed in homogenates",hom_miss)

# Keep the clusters occurring in more than 100 samples
x <- sort(table(D$cluster),decreasing=TRUE)
Fs <- D[D$cluster %in% names(x[x>100]),]
Fs$cluster <- factor(Fs$cluster,levels=names(x[x>100]))

lys_sample_mean <- sum(D$lysate_reads)/length(unique(D$sample))
hom_sample_mean <- sum(D$homogenate_reads)/length(unique(D$sample))

L <- aggregate(Fs$lysate_reads[Fs$lysate_reads!=0]~Fs$cluster[Fs$lysate_reads!=0],FUN=mean)
H <- aggregate(Fs$homogenate_reads[Fs$homogenate_reads!=0]~Fs$cluster[Fs$homogenate_reads!=0],FUN=mean)
#L <- aggregate(F$lysate_reads~F$cluster,FUN=mean)
#H <- aggregate(F$homogenate_reads~F$cluster,FUN=mean)
colnames(L) <- c("cluster","mean_lys_reads")
colnames(H) <- c("cluster","mean_hom_reads")
LH <- merge(L, H, by="cluster")

set1 <- LH$cluster[order(LH$mean_lys_reads/lys_sample_mean+LH$mean_hom_reads/hom_sample_mean, decreasing=TRUE)]
set2 <- LH$cluster[order(LH$mean_lys_reads/LH$mean_hom_reads,decreasing=TRUE)]
set3 <- LH$cluster[order(LH$mean_lys_reads/LH$mean_hom_reads)]
set4 <- LH$cluster[order(LH$mean_lys_reads/lys_sample_mean+LH$mean_hom_reads/hom_sample_mean)]

cutoff <- 50
set1 <- set1[1:cutoff]
set2 <- set2[1:cutoff]
set3 <- set3[1:cutoff]
set4 <- set4[1:cutoff]

cat("\nThe",cutoff,"clusters with the largest mean reads (lysates + homogenates)\n",file=file,append=TRUE)
cat(as.character(set1),sep="\n",file=file,append=TRUE)
tax_stats("clusters with the largest mean reads", set1)

cat("\nThe",cutoff,"clusters with the largest surplus of lysate reads (difference in mean reads)\n",file=file,append=TRUE)
cat(as.character(set2),sep="\n",file=file,append=TRUE)
tax_stats("clusters with the largest lysate surplus", set2)

cat("\nThe",cutoff,"clusters with the largest surplus of homogenate reads (difference in mean reads)\n",file=file,append=TRUE)
cat(as.character(set3),sep="\n",file=file,append=TRUE)
tax_stats("clusters with the largest homogenate surplus", set3)

cat("\nThe",cutoff,"clusters with the smallest mean reads (lysates + homogenates)\n",file=file,append=TRUE)
cat(as.character(set4),sep="\n",file=file,append=TRUE)
tax_stats("clusters with the smallest mean reads", set4)

