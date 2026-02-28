# Symbio-COI-Barcoding
This repository contains scripts and other files used for analysing amplicon sequencing data (COI) for the publication by XXX *et al.* titled *"Accuracy of occurrence and abundance estimates from insect metabarcoding."*

The starting data comprised merged reads from Illumina amplicon sequencing. 

The process is divided into two main steps, separated into the scripts LSD1.PY and LSD2.py. 

## LSD1.py

To run this script, you have to be in the folder with your raw data (fastq.gz files). 

This script uses VSEARCH v 2.22.1 to:
1. Filter fastq.gz files by length (minimum: 308, maximum: 314).
2. Transform files from fastq.gz to fasta format.
3. Dereplicates sequences (merging full-length identical sequences and writing a fasta file with representative sequences for each group, indicating the size of the group in the header).
4. Sorts sequences by group size and removes singletons.
5. Denoises: identifies biological sequences (ZOTUs) by removing sequencing errors and chimeras.
6. Creates a zOTU table.
7. Adds sequences to zOTU table using the custom script *add_seq_to_zotu*, also found on this repository.

All these steps are performed for each sample individually. 

The most relevant outputs for this step are:
- A zOTU fasta file for each sample (with the name pattern: sample.zotus.fasta)
- A zOTU table with sequences for each sample (sample_zotu_table_with_seq.txt)


## LSD2.py

To run this script, you have to provide the path to the folder where your raw data is located. 

The script uses VSEARCH v 2.22.1,  USEARCH v 11.0.667, and custom Python code to perform:
1. Merging of zOTU information from individual samples into a single zOTU table. 
2. OTU picking and chimeric sequence detection.
3. Taxonomic annotation to zOTUs and OTUs.
4. Creation of zOTU and OTU tables with sequences, taxonomic assignment, and abundances.

The most relevant outputs for this step are:
- zOTU reassignment fasta file (from combining data from all samples, new_zotus.fasta)
- Files with taxonomy information (zotus.tax and otus.tax)
- An OTU table (OTU_Table.txt)
- A zOTU table (zotu_table_expanded.txt)


