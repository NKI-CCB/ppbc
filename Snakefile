configfile: "config.yaml"

rule all:
    input:
        "reports/02_QC_salmon.html"
        #expand("data/RNA-seq/salmon/{sample}/quant.sf", sample=config['samples'])

rule fastqc:
    input:
        "data/RAW/{sample}.fastq.gz"
    output:
        html = "results/fastqc/{sample}_fastqc.html",
        zip = "results/fastqc/{sample}_fastqc.zip"
    threads: 5
    conda:
        "environment.yml"         
    shell:
        """
        ~/tools/FastQC/fastqc {input} -o results/fastqc/ -t {threads}
        """

rule salmon_index:
    input:
        "src/01-build-salmon-index.sh"
    output:
        #directory("data/external/index/grch38_index")
        "data/external/index/grch38_index/hash.bin"
    conda:
        "environment.yml"        
    shell:
        "bash {input}"
        
rule salmon_quant:
    input:
        "src/01-salmon.sh"
    output:
        "data/RNA-seq/salmon/{sample}/quant.sf"
    params:
        index = "data/external/index/grch38_index",
        outdir = "data/RNA-seq/salmon/{sample}"
    conda:
        "environment.yml"     
    shell:
        """
        bash {input}
        """

rule multiqc:
    input:
        fastqc = "results/fastqc/",
        salmon_quant = "data/RNA-seq/salmon/"
    output:
        "reports/multiqc_report.html"
    conda:
        "environment.yml"             
    shell:
        """multiqc {input.fastqc} {input.salmon_quant} -o results -n multiqc_report.html --force"""

rule process_metadata:
    input:
        raw_metadata="data/external/Hercoderingslijst_v09032020_KM.xlsx"
    output:
        report="reports/01_process_metadata_tximport.html",
        multiple_patient_fastqs="data/metadata/01_patients_with_multiple_fastqs.csv",
        excluded_samples="data/metadata/01_pre_excluded_samples.csv",
        sample_annot="data/metadata/01_sample_annot.tsv",
        gene_annot="data/metadata/01_tx_annot.tsv",
        tx="data/Rds/01_tx.Rds",
        rdata="reports/01_process_metadata_tximport.RData"
    script:
        "reports/01_process_metadata_tximport.Rmd"
        
rule QC_salmon:
    input:
      #A tx import object
      tx="data/Rds/01_tx.Rds", 
      #Sample and gene annotation
      sample_annot="data/metadata/01_sample_annot.tsv",
      refseq_db="data/external/refseqid_genename_hg38.txt",
      recent_samples="data/metadata/new_samples_jun2019.txt",
      #Two multiqc report files
      alignstats="reports/multiqc_data/multiqc_general_stats.txt",
      or_summary="reports/multiqc_data/mqc_fastqc_overrepresented_sequencesi_plot_1.txt",      
      #Blast files created via the web browser based fastas generated in the Rmd
      blast_or="results/fastqc/overrepresented-BLAST-HitTable.csv",
      blast_failed="results/fastqc/failedor-BLAST-HitTable.csv",
      gc_blast="results/fastqc/gc_or_Alignment-HitTable.csv"
    output:
      "reports/02_QC_salmon.html",
      #Overrepresented fasta sequences aggregated in the report
      "results/fastqc/overrepresented.fa",
      "results/fastqc/failed_overrepresented.fa",
      "results/fastqc/gc_or.fa",
      #Upset plots
      "data/metadata/02_QC_samples_discarded_by_reason.pdf",
      "data/metadata/02_QC_patients_discarded_by_reason.pdf",
      #Kept vs discarded metadata
      "data/metadata/02_discarded_samples.csv",
      "data/metadata/02_sample_annot_filtered.csv",
      #Txdata from only kept samples
      "data/Rds/02_tx_clean.Rds",
      #A dds from only kept samples, with replicates collapsed
      "data/Rds/02_QC_dds.Rds",
      "reports/02_QC_salmon.RData"
    script:
      "reports/02_QC_salmon.Rmd"
    
    
    
        