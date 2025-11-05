# paper-metabar-accuracy
Repository with data and code for the ms evaluating the accuracy of occurrence and abundance estimations from metabarcoding data.

## R scripts

### `Accuracy_paper_map.R`
Generates a map of Sweden with marked sampling points, where 856 homogenate-lysate sample pairs were collected. 
It marks with color the six different habitats.

### `plot_abundance_calib.R`
Generates plots illustrating the power of various calibration procedures using spike-ins in reducing the variance in read numbers. All samples with spike-in reads are used, that is, the lysate data are based on all lysate samples, not only the ones that were also homogenized.
Only samples with all biological spike-ins and all synthetic spike-ins are kept.
The generated plots are `Fig_lysate_cals.jpg`, `Fig_homogenate_cals.jpg`and `Fig_calibrations.jpg`. The last plot combines the two others. The panels illustrate raw lysate reads of biological spikeins, lysate reads of biological spikeins calibrated with the other biological spikeins (leave-one-out), the raw homogenate reads, homogenate reads calibrated with synthetic spike-ins, homogenate reads calibrated with biological spikeins, and homogenate reads calibrated using a linear model trying to optimize the combination of information from biological and synthetic spike-ins.

### `match_lysate_homogenate.R`
Generates files with calibrated lysate, homogenate and matched lysate-homogenate data.

The files are as follows (also provided in the github repo)
1. `shared_lysates.tsv`. Data frame with clusters as rows and lysate reads in each of the shared samples as columns.
2. `shared_homogenates.tsv`. Ditto for homogenates. 
3. `matched_lysate_homogenate.tsv`. Melted data frame with samples, lysate reads and homogenate reads for all cases where a cluster was detected either in homogenates or lysates.

All files contain full taxonomic information for the cluster. The last file is > 50 MB and needs to be gzipped in github to avoid warning messages.

### `plot_matches.R`
Generates plots illustrating the distributions of calibrated read counts for lysates and homogenates for different groups of OTU clusters. The generated figures illustrate the following:

Main figures:
1. `Fig2_All_Samples_500_Ratio.jpg`: Lysate-to-homogenate read ratios for 856 treatment pairs. Part of Figure 2 from the paper.
2. `Fig2_treatment_exclusive_clusters.jpg`: Arthropod clusters detected exclusively in either lysates or homogenates, aggregated by family. Second part of Figure 2 from the paper.

Supplementary figures:
1. `Fig_S2_matched_largest_lysate_surplus.jpg`: Top 20 clusters (among those clusters occurring in more than 100 samples) with largest surplus in lysates (top panel) when compared with homogenates (bottom panel).
2. `Fig_S3_matched_largest_homogenate_surplus.jpg`: Top 20 clusters (among those clusters occurring in more than 100 samples) with largest surplus in homogenates (top panel) when compared with homogenates (bottom panel).
3. `Fig_S4_Scatter_all_famillies.jpg`: Scatterplot summarizing lysate to homogenate number of reads obtained per family (limited to Arthropoda), across 856 samples (treatment pairs).
4. `Fig_S5_15samples_reads_biomass.jpg`: Selected samples - sequencing depth and biomass. Histogram shows the total number of reads per sample for both lysate and homogenate as well as the biomass (wet weight) of each sample.
5. `Fig_S6_15Samples_Top10_Ratio.jpg`: Lysate-to-homogenate read ratios for 15 samples.

And additional exploratory figures:
`Fig_matched_most_occurring.jpg`: Illustrates the pattern for the 20 clusters detected in the largest number of samples among the samples with both lysate and homogenate data.
`Fig_matched_most_reads.jpg`: Pattern for the 20 clusters with the largest read counts among those clusters occurring in more than 100 samples
`Fig_matched_least_reads.jpg`: Pattern for the 20 clusters with the smallest read counts among those clusters occurring in more than 100 samples

### `generate_match_stats.R`
Generates more detailed statistics for the matching of lysate and homogenate reads. The statistics are printed to the file `match_stats.txt`.