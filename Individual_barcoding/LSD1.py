#! /usr/bin/env python3

##filter
for file in *.fastq.gz; do
    SampleName=`basename $file .fastq.gz`
    vsearch --fastx_filter $SampleName.fastq.gz --fastq_minlen 308 --fastq_maxlen 314 --fastqout $SampleName_filtered.fastq.gz
done

##Fastq to fasta
for file in *.fastq; do
    SampleName=`basename $file .fastq`
    vsearch -fastq_filter $SampleName.fastq -fastaout $SampleName.fasta -relabel "$SampleName"._ -fasta_width 0
done

##dereplicating
for file in *.fasta; do
    SampleName=`basename $file .fasta`
    vsearch -derep_fulllength $SampleName.fasta -output "$SampleName".derep.fasta -sizeout -uc "$SampleName".derep_info.txt
done

##picking representative sequence
for file in *derep.fasta; do
    SampleName=`basename $file .derep.fasta`
    vsearch -sortbysize $SampleName.derep.fasta --output $SampleName.sorted.fasta -minsize 2
done

mkdir derep && mv *derep* derep/

##denoising
for file in *sorted.fasta; do
        SampleName=`basename $file .sorted.fasta`
        usearch -unoise3 $SampleName.sorted.fasta -zotus $SampleName.zotus.fasta -tabbedout $SampleName.denoising.summary.txt -minsize 1
done

mkdir sorted && mv *sorted.fasta sorted/
mkdir denoising_summary && mv *denoising.summary.txt denoising_summary/

##creating zotu tables
for file in *.fasta; do
    SampleName=`basename $file .zotus.fasta`
    usearch -otutab ./filtered/"$SampleName".filtered.fasta -zotus $SampleName.zotus.fasta -otutabout "$SampleName"_zotu_table.txt -threads 60
done

##Add sequence to zOTU table
for file in *.fasta; do
    SampleName=`basename $file .zotus.fasta`
    python3 /home/karol.nowak/ela_iwaszk/add_seq_to_zotu.py "$SampleName"_zotu_table.txt "$SampleName".zotus.fasta "$SampleName"_zotu_table_with_seq.txt
done



