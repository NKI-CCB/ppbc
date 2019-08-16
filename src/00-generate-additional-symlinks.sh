#!/bin/bash

#for f in 'rHLE-FA-156.190621.HiSeq4000.FCB.lane5.gcap_19_01.R1.fastq.gz' 'rHLE-FA-286.190621.HiSeq4000.FCB.lane6.gcap_19_01.R1.fastq.gz' 'rHLE-FA-163.190621.HiSeq4000.FCB.lane6.gcap_19_01.R1.fastq.gz' 'rHLE-FA-308.190621.HiSeq4000.FCB.lane5.gcap_19_01.R1.fastq.gz'; do echo "/DATA/share/postpartbc/fastq_jun2019/$f"; done

for f in 'rHLE-FA-156.190621.HiSeq4000.FCB.lane5.gcap_19_01.R1.fastq.gz' 'rHLE-FA-286.190621.HiSeq4000.FCB.lane6.gcap_19_01.R1.fastq.gz' 'rHLE-FA-163.190621.HiSeq4000.FCB.lane6.gcap_19_01.R1.fastq.gz' 'rHLE-FA-308.190621.HiSeq4000.FCB.lane5.gcap_19_01.R1.fastq.gz'; do ln -s "/DATA/share/postpartbc/fastq_jun2019/$f" "/DATA/share/postpartumbc/data/RAW/$f"; done
