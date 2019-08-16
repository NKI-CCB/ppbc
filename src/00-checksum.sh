#!/bin/bash

#for file in /DATA/share/postpartbc/fastq_jun2019/*/runs/*/v*/result/*R1.fastq.gz; do echo $file; done

for file in /DATA/share/postpartbc/fastq_jun2019/*/runs/*/v*/result/*R1.fastq.gz; do echo "$file"; md5sum $file >> /DATA/share/postpartumbc/data/metadata/md5_jun2019.txt; done
