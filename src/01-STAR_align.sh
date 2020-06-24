#for f in `ls /DATA/share/postpartumbc/data/RAW/fastq`; do echo ${f}; done
#for f in `ls /DATA/share/postpartumbc/data/RAW/fastq`; do prefix="${f%.R1.fastq.gz}"; echo "${prefix}"; done


#cd /DATA/share/postpartumbc/data/STAR/run2

for f in `ls /DATA/share/postpartumbc/data/RAW/*.R1.fastq.gz`
do
    basefile="$(basename -- $f)"
    prefix="${basefile%.R1.fastq.gz}"
    echo "Starting ${prefix}..."
STAR --runThreadN 16 \
  --genomeDir /DATA/share/postpartumbc/data/external/index/STAR_grch38_index/ \
  --readFilesIn "${f}" \
  --outFileNamePrefix "/DATA/share/postpartumbc/data/RNA-seq/star/${prefix}" \
  --outFilterMultimapNmax 20 \
  --outFilterMismatchNmax 1 \
  --outSAMmultNmax 1 \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts \
  --sjdbGTFfile /DATA/share/postpartumbc/data/external/index/Homo_sapiens.GRCh38.94.gtf \
  --readFilesCommand zcat \
  --outWigType None
done
