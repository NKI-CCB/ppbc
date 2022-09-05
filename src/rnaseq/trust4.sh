#for f in `ls $PWD/data/rnaseq/fastq`; do echo ${f}; done
#for f in `ls $PROJDIR/data/rnaseq/fastq/*.R1.fastq.gz`; do basefile="$(basename -- $f)"; prefix="${basefile%%.*}"; echo "Starting ${prefix}"; done

#Usage: ./run-trust4 [OPTIONS]
#	Required:
#		-b STRING: path to bam file
#		-1 STRING -2 STRING: path to paired-end read files
#		-u STRING: path to single-end read file
#		-f STRING: path to the fasta file coordinate and sequence of V/D/J/C genes
#	Optional:
#		--ref STRING: path to detailed V/D/J/C gene reference file, such as from IMGT database. (default: not used). (recommended) 
#		-o STRING: prefix of output files. (default: inferred from file prefix)
#		-t INT: number of threads (default: 1)
#		--barcode STRING: if -b, bam field for barcode; if -1 -2/-u, file containing barcodes (defaul: not used)
#		--barcodeRange INT INT CHAR: start, end(-1 for lenght-1), strand in a barcode is the true barcode (default: 0 -1 +)
#		--abnormalUnmapFlag: the flag in BAM for the unmapped read-pair is nonconcordant (default: not set)
#		--noExtraction: directly use the files from provided -1 -2/-u to assemble (default: extraction first)
#		--stage INT: start TRUST4 on specified stage (default: 0)
#			0: start from beginning (candidate read extraction)
#			1: start from assembly
#			2: start from annotation
#			3: start from generating the report table

#Example
#./run-trust4 -b example/example.bam -f hg38_bcrtcr.fa --ref human_IMGT+C.fa

#Test output
#f=$PWD/data/rnaseq/fastq/rHLE-FA-100.180816.HiSeq4000.FCB.lane5.gcap_17_11.R1.fastq.gz

PROJDIR=$PWD
OUTDIR=$PROJDIR/data/TRUST
TRUSTDIR=$PROJDIR/bin/TRUST4

mkdir -p $OUTDIR

for f in `ls $PROJDIR/data/rnaseq/fastq/*.R1.fastq.gz`
do
    basefile="$(basename -- $f)"
    prefix="${basefile%%.*}"
    
    mkdir -p $OUTDIR
    cd $OUTDIR
    
    echo "Starting ${prefix}..."
    $TRUSTDIR/run-trust4 -u $f -t 16 -f "$TRUSTDIR/hg38_bcrtcr.fa" --ref "$TRUSTDIR/human_IMGT+C.fa" -o $prefix
      
done
