#!/bin/bash

#for file in /DATA/share/postpartumbc/data/RAW/fastq_jun2019/*fastq.gz; do salmon quant -i data/external/index/grch38_index -l A --gcBias --seqBias --posBias -r $file -o '{params.outdir}'; done
#for file in /DATA/share/postpartumbc/data/RAW/fastq_jun2019/*fastq.gz; do samp=`basename ${file}`; outname=${samp%.fastq*}; echo "Processing sample ${outname}";  done

for file in /DATA/share/postpartumbc/data/RAW/fastq_jun2019/*fastq.gz;
do 
	samp=`basename ${file}`
	outname=${samp%.fastq*}
	echo "Processing sample ${outname}"
	salmon quant -i /DATA/share/postpartumbc/data/external/index/grch38_index -l A --gcBias --validateMappings --seqBias --posBias -p 32 -r $file -o "/DATA/share/postpartumbc/data/RNA-seq/salmon/${outname}"

done

#The version of salmon has changed, so we also realign the older files
for file in /DATA/share/postpartumbc/data/RAW/fastq/*fastq.gz;
do 
	samp=`basename ${file}`
	outname=${samp%.fastq*}
	echo "Processing sample ${outname}"
	salmon quant -i /DATA/share/postpartumbc/data/external/index/grch38_index -l A --gcBias --seqBias --validateMappings --posBias -p 32 -r $file -o "/DATA/share/postpartumbc/data/RNA-seq/salmon/${outname}"

done
