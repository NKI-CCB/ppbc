configfile: "config.yaml"

#### Conda and Snakemake ####
#To activate conda environment, use: conda env create -f envs/environment.yml
#If environment.yml has been updated, update the environment with: conda env update --prefix ./envs --file environment.yml  --prune
#Save environment to text with: conda env export > envs/environment.yml
#For condensed version: conda env export --from-history
#To use a conda envinroment with snakemake: snakemake -n --use-conda
#To create the environments without running any rules: snakemake -n --use-conda --create-envs-only

rule all:
  input:
    "reports/06_diffex_lrt.html"
    #"dag.svg"
    #expand("data/RNA-seq/salmon/{sample}/quant.sf", sample=config['samples'])

rule fastqc:
  input:
    expand("data/RAW/{sample}.fastq.gz", sample=config['samples'])
  output:
    html = expand("results/fastqc/{sample}_fastqc.html", sample=config['samples']),
    zip = expand("results/fastqc/{sample}_fastqc.zip", sample=config['samples'])
  threads: 5
  conda:
    "envs/environment.yml"         
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
    "envs/environment.yml"        
  shell:
    "bash {input}"

rule salmon_quant:
  input:
    "src/01-salmon.sh"
  output:
    expand("data/RNA-seq/salmon/{sample}/quant.sf", sample=config['samples'])
  params:
    index = "data/external/index/grch38_index",
    outdir = expand("data/RNA-seq/salmon/{sample}", sample=config['samples'])
  conda:
    "envs/environment.yml"     
  shell:
    """
    bash {input}
    """

rule multiqc:
  input:
    fastqc = "results/fastqc/",
    salmon_quant = "data/RNA-seq/salmon/"
  output:
    report="reports/multiqc_report.html",
    #alignstats="reports/multiqc_data/multiqc_general_stats.txt", #Decouple
    #or_summary="reports/multiqc_data/mqc_fastqc_overrepresented_sequencesi_plot_1.txt" #Decouple
  conda:
    "envs/environment.yml"             
  shell:
    """multiqc {input.fastqc} {input.salmon_quant} -o results -n multiqc_report.html --force"""

rule process_metadata:
  input:
    #expand("data/RNA-seq/salmon/{sample}/quant.sf", sample=config['samples']), #Decouple
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

rule color_palettes:
  input:
    rmd="reports/color_palettes.Rmd",
    script="src/rmarkdown.R",
    dds="data/Rds/03_dds_PAM50.Rds"
  output:
    html="reports/color_palettes.html",
    color_palettes="data/Rds/color_palettes.Rds"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

rule survival_colors:
  output:
    "data/Rds/survival_colors.Rds"
  shell:
    "Rscript src/survival_colors.R"

rule surv_est:
  input:
    script="src/rmarkdown.R",
    rmd="reports/04_survival_and_ESTIMATE.Rmd",
    dds="data/Rds/03_dds_PAM50.Rds",
    tx_annot="data/metadata/01_tx_annot.tsv",
    color_palette="data/Rds/color_palettes.Rds"
  output:
    #Files from ESTIMATE
    "data/RNA-seq/hugo_fpkm.txt",
    "results/ESTIMATE/filterCommonGenes.gct",
    "results/ESTIMATE/results_estimate_score.gct",
    #Survival metadata
    "data/metadata/04_samples_with_missing_survival_data.csv",
    "data/metadata/04_survival_metadata.csv",
    "data/metadata/04_samples_excluded_survival.xlsx",
    #Samples x features matrix for Cox regressions
    "data/Rds/04_survdata.Rds",
    #Kaplan-meier survival curves
    "results/survival/04_kaplan_meiers.pdf",
    #Updated colData for dds
    "data/metadata/04_sample_annot_filtered_PAM50_EST.csv",
    #Updated dds
    "data/Rds/04_dds_PAM50_EST.Rds",
    #Notebook environment
    "reports/04_survival_and_ESTIMATE.RData",
    html="reports/04_survival_and_ESTIMATE.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

pca_pdfs = [
  "scree", "sigPCA_batch", "firstPCA_batch",
  "sigPCA_PPBC", "firstPCA_PPBC", "sigPCA_Pam50"
]
    
rule batch_effects:
  input:
    script="src/rmarkdown.R",
    rmd="reports/05_batch_effects.Rmd",
    dds="data/Rds/04_dds_PAM50_EST.Rds",
    gx_annot="data/metadata/01_tx_annot.tsv",
    cp="data/Rds/color_palettes.Rds"
  output:
    expand("results/PCA/05_{pca}.pdf", pca=pca_pdfs),
    metadata="data/metadata/05_sample_annot_filtered.csv",
    dds="data/Rds/05_dds_PAM50_batch.Rds",
    rdata="reports/05_batch_effects.RData",
    html="reports/05_batch_effects.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

#Setting the number of threads prevents unwanted parallelization    
rule lrt_diffex:
  input:
    dds="data/Rds/05_dds_PAM50_batch.Rds"
  output:
    "data/Rds/06_vsd.Rds",
    "data/Rds/06_vsd_nolac.Rds",
    "data/Rds/06_ddsLRT.Rds",
    "data/Rds/06_ddsLRT_nolac.Rds",
    "data/Rds/06_ddsLRT_pam.Rds",
    "data/Rds/06_ddsLRT_pam_nolac.Rds",
    "data/Rds/06_ddsLRT_batch.Rds",
    "data/Rds/06_ddsLRT_batch_nolac.Rds"
  shell:
    """
    export OMP_NUM_THREADS=1 
    Rscript src/06_diffex_lrt.R
    """

dds_lrt = [
    "vsd", "vsd_nolac", "ddsLRT",
    "ddsLRT_nolac", "ddsLRT_pam", "ddsLRT_pam_nolac",
    "ddsLRT_batch", "ddsLRT_batch_nolac"
    ]
    
gene_sets = [
  "c5.bp.v7.0.symbols", "h.all.v7.0.symbols",
  "c2.cp.v7.0.symbols", "c2.cgp.v7.0.symbols"
]

lrt_reports = [
  "LRT", "LRT.nolac", "Batch",
  "Batch.nolac", "Pam", "Pam.nolac"
]

rule lrt_report:
  input:
    expand("data/Rds/06_{dds}.Rds", dds=dds_lrt),
    expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    dds="data/Rds/05_dds_PAM50_batch.Rds",
    #ImmPort database: https://www.innatedb.com/redirect.do?go=resourcesGeneLists
    immune_genes="data/external/gene-sets/InnateDB_genes.csv",
    gx_annot="data/metadata/01_tx_annot.tsv",
    cp="data/Rds/color_palettes.Rds",
    tools="src/deseq_report_functions.R",
    rmd="reports/06_diffex_lrt.Rmd",
    script="src/rmarkdown.R"
  output:
    padj_comp="results/diffex/figs/06_LRT/06_allgenes_fdrvsexcludelac.pdf",
    padj_comp_sig="results/diffex/figs/06_LRT/06_siggenes_fdrvsexcludelac.pdf",
    sig_genes="results/diffex/06_LRT_sig_genes.xlsx",
    all_genes="results/diffex/06_LRT_allgenes.xlsx",
    pathways="results/diffex/06_LRT_pathways.xlsx",
    heatmaps=expand("results/diffex/figs/06_LRT/heatmaps/hm{comp}.pdf", comp=lrt_reports),
    volcanos=expand("results/diffex/figs/06_LRT/volcano_plots/volc{comp}.outliers.jpeg", comp=lrt_reports),
    #rdata="reports/06_diffex_lrt.Rdata",
    html="reports/06_diffex_lrt.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

pairwise_refs = ["non_prbc", "ppbc_lac", "prbc"]

pairwise_apeglm = [
  "prbc_vs_non_prbc", "ppbc_lac_vs_non_prbc", "ppbc_inv_vs_non_prbc", 
  "prbc_vs_ppbc_lac", "ppbc_inv_vs_prbc", "ppbc_inv_vs_ppbc_lac"
]
    
rule diffex_pairwise:
  input:
    dds="data/Rds/05_dds_PAM50_batch.Rds"
  output:
    pairwise_dds=expand("data/Rds/07_dds_pairwise_ref_{pw}.Rds", pw=pairwise_refs),
    apeglm_results=expand("data/Rds/07_ape_{pw}.Rds", pw=pairwise_apeglm)
  shell:
    """
    export OMP_NUM_THREADS=1 
    Rscript src/07_diffex_pairwise.R
    """
    
pairwise_prefixes = [
  "inv_lac", "inv_nonprbc", "inv_prbc",
  "lac_nonprbc", "lac_prbc", "prbc_nonprbc"
]

rule pairwise_report:
  input:
    pairwise_dds=expand("data/Rds/07_dds_pairwise_ref_{pw}.Rds", pw=pairwise_refs),
    sets=expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    immune_genes="data/external/gene-sets/InnateDB_genes.csv",
    cp="data/Rds/color_palettes.Rds",
    sp="data/Rds/survival_colors.Rds",
    vsd="data/Rds/06_vsd.Rds",
    tools="src/deseq_report_functions.R",
    apeglm_results=expand("data/Rds/07_ape_{pw}.Rds", pw=pairwise_apeglm),
    gx_annot="data/metadata/01_tx_annot.tsv",
    rmd="reports/07_diffex_pairwise.Rmd",
    script="src/rmarkdown.R"
  output:
    "results/diffex/figs/07_pairwise/heatmaps/n_siggenes_pairwise_hm.pdf",
    #"results/diffex/figs/07_pairwise/heatmaps/07_commoninvpaths_hm.pdf",
    "results/diffex/figs/07_pairwise/upset_pairwise.pdf",
    "results/diffex/07_pairwise_comparisons_allgenes.xlsx",
    "results/diffex/07_pairwise_comparisons_sig_genes.xlsx",
    "results/diffex/07_pairwise_comparisons_pathways.xlsx",
    expand("results/diffex/figs/07_pairwise/heatmaps/hm_{pf}.pdf", pf=pairwise_prefixes),
    expand("results/diffex/figs/07_pairwise/volcano_plots/volc_{pf}.jpeg", pf=pairwise_prefixes),
    expand("results/diffex/figs/07_pairwise/volcano_plots/volc_{pf}.outliers.jpeg", pf=pairwise_prefixes),
    html="reports/06_diffex_pairwise.html"
  shell:
     "Rscript {input.script} {input.rmd} $PWD/{output.html}"

ovr_comps = ["inv_vs_rest", "prbc_vs_rest", "lac_vs_rest", "nonprbc_vs_rest"]
     
rule diffex_one_vs_rest:
  input:
    dds="data/Rds/05_dds_PAM50_batch.Rds"
  output:
    dds_ovr=expand("data/Rds/08_dds_ovr_{comp}.Rds", comp=ovr_comps),
    ape_ovr=expand("data/Rds/08_ape_ovr_{comp}.Rds", comp=ovr_comps)
  shell:
    """
    export OMP_NUM_THREADS=1 
    Rscript src/08_diffex_onevsrest.R
    """
    
rule workflow_diagram:
  conda:
    "envs/environment.yml"
  shell:
    "snakemake --use-conda --dag | dot -Tsvg > dag.svg"
    #To view via shell: display dag.svg
    
