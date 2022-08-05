#!/bin/bash


#fastqc /DATA/share/postpartumbc/data/RAW/*.fastq.gz -o /DATA/share/postpartumbc/results/fastqc

#Additional files that wree delivered later
#for f in 'rHLE-FA-156.190621.HiSeq4000.FCB.lane5.gcap_19_01.R1.fastq.gz' 'rHLE-FA-286.190621.HiSeq4000.FCB.lane6.gcap_19_01.R1.fastq.gz' 'rHLE-FA-163.190621.HiSeq4000.FCB.lane6.gcap_19_01.R1.fastq.gz' 'rHLE-FA-308.190621.HiSeq4000.FCB.lane5.gcap_19_01.R1.fastq.gz'; do fastqc /DATA/share/postpartumbc/data/RAW/$f -o /DATA/share/postpartumbc/results/fastqc; done

multiqc /DATA/share/postpartumbc/results/fastqc/ /DATA/share/postpartumbc/data/RNA-seq/salmon -o /DATA/share/postpartumbc/reports


