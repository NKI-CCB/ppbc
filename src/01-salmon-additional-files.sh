#!/bin/bash

#A number of samples had to be resequenced and were delivered later than the rest

cd /DATA/share/postpartbc/fastq_jun2019

for file in 'rHLE-FA-156.190621.HiSeq4000.FCB.lane5.gcap_19_01.R1.fastq.gz' 'rHLE-FA-286.190621.HiSeq4000.FCB.lane6.gcap_19_01.R1.fastq.gz' 'rHLE-FA-163.190621.HiSeq4000.FCB.lane6.gcap_19_01.R1.fastq.gz' 'rHLE-FA-308.190621.HiSeq4000.FCB.lane5.gcap_19_01.R1.fastq.gz';
do 
	samp=`basename ${file}`
	outname=${samp%.fastq*}
	echo "Processing sample ${outname}"
	salmon quant -i /DATA/share/postpartumbc/data/external/index/grch38_index -l A --gcBias --validateMappings --seqBias --posBias -p 32 -r $file -o "/DATA/share/postpartumbc/data/RNA-seq/salmon/${outname}"

done
