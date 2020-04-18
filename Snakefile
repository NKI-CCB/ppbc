configfile: "config.yaml"

rule all:
  input:
    "reports/04_survival_and_ESTIMATE.html"
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
        fastqc {input} -o results/fastqc/ -t {threads}
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

#Rmarkdown files may only have a single output file in Snakemake.
#We'll keep track of other files, but only the html report can be used as a rule target.

rule process_metadata:
  input:
    raw_metadata="data/external/Hercoderingslijst_v09032020_KM.xlsx",
    rmd="reports/01_process_metadata_tximport.Rmd",
    script="src/rmarkdown.R"
  output:
    html="reports/01_process_metadata_tximport.html",
    multiple_patient_fastqs="data/metadata/01_patients_with_multiple_fastqs.csv",
    excluded_samples="data/metadata/01_pre_excluded_samples.csv",
    sample_annot="data/metadata/01_sample_annot.tsv",
    gene_annot="data/metadata/01_tx_annot.tsv",
    tx="data/Rds/01_tx.Rds",
    rdata="reports/01_process_metadata_tximport.RData"
  shell:
    #"mkdir -p reports\n"
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
  

rule QC_salmon:
  input:
    script="src/rmarkdown.R",
    rmd="reports/02_QC_salmon.Rmd",
    #A tx import object
    tx="data/Rds/01_tx.Rds", 
    #Sample and gene annotation
    sample_annot="data/metadata/01_sample_annot.tsv",
    refseq_db="data/external/refseqid_genename_hg38.txt",
    recent_samples="data/metadata/new_samples_jun2019.txt",
    #Two multiqc report files
    alignstats="reports/multiqc_data/multiqc_general_stats.txt",
    or_summary="reports/multiqc_data/mqc_fastqc_overrepresented_sequencesi_plot_1.txt",      
    #Blast files created via the web browser based on fastas generated in the Rmd
    blast_or="results/fastqc/overrepresented-BLAST-HitTable.csv",
    blast_failed="results/fastqc/failedor-BLAST-HitTable.csv",
    gc_blast="results/fastqc/gc_or_Alignment-HitTable.csv"
  output:
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
    "reports/02_QC_salmon.RData",
    html="reports/02_QC_salmon.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

rule pam50_ihc:
  input:
    script="src/rmarkdown.R",
    rmd="reports/03_PAM50_IHC.Rmd",
    #A DESeqDataSet
    dds="data/Rds/02_QC_dds.Rds"
  output:
    "data/metadata/03_ihc_outliers.csv",
    # Samples discarded due to clear incongruity
    "data/metadata/03_removed_pam50_outliers.csv",
    # New dds post Pam50 determination
    "data/Rds/03_dds_PAM50.Rds",
    # New metadata post PAM50 determination
    "data/metadata/03_sample_annot_filtered_PAM50.csv",
    html="reports/03_PAM50_IHC.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

kaplan_meiers = [
  "drs_breastfeeding", "drs_months_inv", "drs_ppbc",
  "km_breastfeeding", "km_months_inv", "km_ppbc",
  "drs_immune_score", "drs_pam50", "drs_stromal_score",
  "km_immune_score", "km_pam50", "km_stromal_score"
  ]

rule surv_est:
  input:
    script="src/rmarkdown.R",
    rmd="reports/04_survival_and_ESTIMATE.Rmd",
    dds="data/Rds/03_dds_PAM50.Rds",
    tx_annot="data/metadata/01_tx_annot.tsv"
  output:
    #Files from ESTIMATE
    "data/RNA-seq/hugo_fpkm.txt",
    "results/ESTIMATE/filterCommonGenes.gct"
    "results/ESTIMATE/results_estimate_score.gct"
    #Survival metadata
    "data/metadata/04_samples_with_missing_survival_data.csv"
    "data/metadata/04_survival_metadata.csv",
    "data/metadata/04_samples_excluded_survival.xlsx",
    #Samples x features matrix for Cox regressions
    "data/Rds/04_survdata.Rds",
    expand("results/survival/04_{km}.pdf", km=kaplan_meiers),
    #"results/survival/04_km_ppbc.pdf"
    #"results/survival/04_km_pam50.pdf"
    #"results/survival/04_km_immune_score.pdf"
    #"results/survival/04_km_stromal_score.pdf"
    #"results/survival/04_drs_ppbc.pdf"
    #"results/survival/04_drs_immune_score.pdf"
    #"results/survival/04_drs_stromal_score.pdf"
    #"results/survival/04_km_months_inv.pdf"
    #"results/survival/04_km_breastfeeding.pdf"
    #"results/survival/04_drs_months_inv.pdf"
    #"results/survival/04_drs_breastfeeding.pdf"
    html="reports/04_survival_and_ESTIMATE.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
