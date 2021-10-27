import pandas as pd
from pathlib import Path

configfile: "config.yaml"

# The Snakefile is given a .py extension so that we can see syntax highlighting in IDEs

# general use:
# -j designates thread number
# -s indicates non-default Snakefile name
# Use -np for dry run
# snakemake {rule or file output} -s Snakefile.py -j 8
# Example:
# snakemake -np gene_reports -s Snakefile.py -j 8

#### Conda and Snakemake ####
#To activate conda environment, use: conda env create -f envs/environment.yml
#If environment.yml has been updated, update the environment with: conda env update --prefix ./envs --file environment.yml  --prune
#Save environment to text with: conda env export > envs/environment.yml
#For condensed version: conda env export --from-history
#To use a conda envinroment with snakemake: snakemake -n --use-conda
#To create the environments without running any rules: snakemake -n --use-conda --create-envs-only

rule all:
  input:
    "reports/16_gene_unity_setup.html"
    #expand("data/RNA-seq/salmon/{sample}/quant.sf", sample=config['samples'])

### RNAseq analyses ###

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

rule star_index:
  input:
    script="src/01-build-STAR-index.sh"
  output:
    directory("data/external/index/STAR_grch38_index")
  params:
    genomeDir="data/external/index/STAR_grch38_index/",
    genomeFastaFiles="data/external/index/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    sjdbGTFfile="data/external/index/Homo_sapiens.GRCh38.94.gtf",
    sjdbOverhang=50 #Read length is 51, ideal value is length - 1
  conda:
    "envs/environment.yml"
  shell:
    """
    mkdir {params.genomeDir}
    bash {input.script}
    """

rule star_align:
  input:
    smpls=expand("data/RAW/{sample}.fastq.gz", sample=config['samples']),
    gtf="data/external/index/Homo_sapiens.GRCh38.94.gtf",
    script="src/01-STAR_align.sh"
  output:
    expand("data/RNA-seq/star/{sample}/Aligned.sortedByCoord.out.bam", sample=config['samples'])
  conda:
    "envs/environment.yml"
  threads: 16
  params:
    genomedir="data/external/index/STAR_grch38_index/",
    outfile=expand("data/RNA-seq/star/{sample}", sample=config['samples']),
    #max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
    multimapnmax=20, #default 10
    #maximum number of mismatches per pair, large number switches off this filter
    mismatchnmax=1, #default 10
    #int: max number of multiple alignments for a read that will be output to the
    #SAM/BAM files. Note that if this value is not equal to -1, the top scoring
    #alignment will be output first
    sammultnmax=1
  shell:
    """
    {input.script}
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
    
rule min_count_threshold:
  input:
    script = "src/rmarkdown.R",
    rmd = "reports/05b_minimum_count_threshold.Rmd",
    dds = "data/Rds/05_dds_PAM50_batch.Rds"
  output:
    dds = "data/Rds/05b_dds_filtered.Rds",
    html = "reports/05b_minimum_count_threshold.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

#Setting the number of threads prevents unwanted parallelization    
rule lrt_diffex:
  input:
    dds="data/Rds/05b_dds_filtered.Rds"
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
    dds="data/Rds/05b_dds_filtered.Rds",
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
    dds="data/Rds/05b_dds_filtered.Rds"
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
    html="reports/07_diffex_pairwise.html"
  shell:
     "Rscript {input.script} {input.rmd} $PWD/{output.html}"

ovr_comps = ["inv_vs_rest", "prbc_vs_rest", "lac_vs_rest", "nonprbc_vs_rest"]
     
rule diffex_one_vs_rest:
  input:
    dds="data/Rds/05b_dds_filtered.Rds"
  output:
    dds_ovr=expand("data/Rds/08_dds_ovr_{comp}.Rds", comp=ovr_comps),
    ape_ovr=expand("data/Rds/08_ape_ovr_{comp}.Rds", comp=ovr_comps)
  shell:
    """
    export OMP_NUM_THREADS=1 
    Rscript src/08_diffex_onevsrest.R
    """
ovr_reports = ["inv_rest", "prbc_rest", "lac_rest", "nonprbc_rest"]

rule one_vs_rest_report:
  input:
    ape_ovr=expand("data/Rds/08_ape_ovr_{comp}.Rds", comp=ovr_comps),
    dds_ovr=expand("data/Rds/08_dds_ovr_{comp}.Rds", comp=ovr_comps),
    sets=expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    immune_genes="data/external/gene-sets/InnateDB_genes.csv",
    cp="data/Rds/color_palettes.Rds",
    sp="data/Rds/survival_colors.Rds",
    vsd="data/Rds/06_vsd.Rds",
    tools="src/deseq_report_functions.R",
    gx_annot="data/metadata/01_tx_annot.tsv",
    rmd="reports/08_diffex_onevsrest.Rmd",
    script="src/rmarkdown.R"
  output:
    "results/diffex/figs/08_onevsrest/upset_onevsrest.pdf",
    "results/diffex/figs/08_onevsrest/08_barplot_onevsrest.jpg",
    "results/diffex/08_one_vs_rest_allgenes.xlsx",
    "results/diffex/08_one_vs_rest_sig_genes.xlsx",
    "results/diffex/08_one_vs_rest_pathways.xlsx",
    expand("results/diffex/figs/08_onevsrest/heatmaps/hm_{ovr}.pdf", ovr=ovr_reports),
    expand("results/diffex/figs/08_onevsrest/volcano_plots/volc_{ovr}.jpeg", ovr=ovr_reports),
    expand("results/diffex/figs/08_onevsrest/volcano_plots/volc_{ovr}.outliers.jpeg", ovr=ovr_reports),
    html="reports/08_diffex_onevsrest.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

diffex_results = [
  "08_one_vs_rest_allgenes",
  "07_pairwise_comparisons_allgenes",
  "06_LRT_allgenes"
  ]

rule genewise_diffex_reports:
  input:
    diffex_results=expand("results/diffex/{result}.xlsx", result=diffex_results),
    cp="data/Rds/color_palettes.Rds",
    sp="data/Rds/survival_colors.Rds",
    tools="src/deseq_report_functions.R",
    vsd="data/Rds/08_vsd_ovr.Rds",
    gx_annot="data/metadata/01_tx_annot.tsv",
    rmd="reports/09_genewise_diffex_reports.Rmd",
    script="src/rmarkdown.R"
  output:
    "results/diffex/figs/milk_vs_IG_genes.pdf",
    html="reports/09_genewise_diffex_reports.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    
rule cibersortX:
  input:
    #Results obtained by uploading hugo_fpkm.txt to the cibersort website with the enumerated parameters
    ciber_res="results/cibersortX/CIBERSORTx_Job3_Results.csv",
    geneEx="data/RNA-seq/hugo_fpkm.txt",
    surv="data/Rds/04_survdata.Rds",
    metadata="data/metadata/05_sample_annot_filtered.csv",
    rmd="reports/10_CibersortX.Rmd",
    script="src/rmarkdown.R"
  output:
    pdf="reports/10_CibersortX.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.pdf}"

#Only those results which include an inv vs something comparison
inv_comps = [
  "07_pairwise_comparisons_allgenes",
  "07_pairwise_comparisons_sig_genes",
  "08_one_vs_rest_allgenes",
  "08_one_vs_rest_sig_genes"
]

rule inv_clustering:
  input:
    expand("results/diffex/{inv_comp}.xlsx", inv_comp=inv_comps),
    rmd="reports/11_clustering_involution.Rmd",
    script="src/rmarkdown.R",
    cp="data/Rds/color_palettes.Rds",
    sp="data/Rds/survival_colors.Rds",
    tools="src/deseq_report_functions.R",
    vsd="data/Rds/08_vsd_ovr.Rds",
    gx_annot="data/metadata/01_tx_annot.tsv",
    surv="data/Rds/04_survdata.Rds"
  output:
    "results/clustering/11_hm_clust_DEG_inv_vs_rest.pdf",
    "results/clustering/11_barplots_ig_clusters.pdf",
    "results/clustering/11_inv_clusters.xlsx",
    "data/Rds/11_ig_clusters.Rds",
    "data/Rds/11_ig_survdata.Rds",
    "results/survival/11_IG_clusters_by_study_group.pdf"
    "results/clustering/11_hm_involution_only_heatmaps.pdf",
    html="reports/11_clustering_involution.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.pdf}"

genewise_cox =[
  "multi_genewise_os", "uni_genewise_os",
  "multi_genewise_drs", "uni_genewise_drs",
  "inv_multi_genewise_os", "inv_uni_genewise_os",
  "inv_multi_genewise_drs", "inv_uni_genewise_drs"
  ]

rule genewise_survival:
  input:
    gx_annot="data/metadata/01_tx_annot.tsv",
    survdata="data/Rds/04_survdata.Rds",
    dds = "data/Rds/08_dds_ovr_inv_vs_rest.Rds",
  output:
    expand("data/Rds/12_{res}.Rds", res=genewise_cox),
    coxdata="data/Rds/12_coxdata.Rds",
    invcoxdata="data/Rds/12_invdata.Rds"
  shell:
    """
    Rscript src/12_genewise_survival.R
    """

rule report_genewise_survival:
  input:
    "src/general_R_tools.R",
    rds=expand("data/Rds/12_{c}.Rds", c=genewise_cox),
    #cox=genewise_cox,
    cp="data/Rds/color_palettes.Rds",
    sp="data/Rds/survival_colors.Rds",
    sets=expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    script="src/12_batch_survival_reports.R",
    rmd="reports/12_genewise_survival.Rmd",
    gx_annot="data/metadata/01_tx_annot.tsv",
    coxdata="data/Rds/12_coxdata.Rds",
    tools="src/enrichment-analysis-functions.R",
    pw="results/diffex/07_pairwise_comparisons_allgenes.xlsx",
    ovr="results/diffex/08_one_vs_rest_allgenes.xlsx"
  output:
    html=expand("reports/12_{rep}.html", rep=genewise_cox),
    csv=expand("results/survival/12_{rep}.csv", rep=genewise_cox)
  shell:
    """
    Rscript {input.script}
    """

rule aggregate_genewise_survival:
  input:
    csv=expand("results/survival/12_{rep}.csv", rep=genewise_cox),
    script="src/rmarkdown.R",
    rmd="reports/12b_aggregate_genewise_survival.Rmd"
  output:
    html="reports/12b_aggregate_genewise_survival.html",
    coxres="results/survival/12_cox_allgenes.xlsx"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

interaction_models = [
  "uni_interaction_os",
  "uni_interaction_drs",
  "multi_interaction_os",
  "multi_interaction_drs"
]    
    
rule surv_inv_int:
  input:
    gx_annot="data/metadata/01_tx_annot.tsv",
    coxdata="data/Rds/12_coxdata.Rds",
    script="src/13_survival_involution_interaction.R"
  output:
    expand("data/Rds/13_{m}.Rds", m=interaction_models)
  shell:
    """
    Rscript {input.script}
    """

rule report_interaction_survival:
  input:
    "src/general_R_tools.R",
    dds="data/Rds/08_dds_ovr_inv_vs_rest.Rds",
    rds=expand("data/Rds/13_{c}.Rds", c=interaction_models),
    cp="data/Rds/color_palettes.Rds",
    sp="data/Rds/survival_colors.Rds",
    sets=expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    script="src/13_batch_involution_interaction_reports.R",
    rmd="reports/13_involutionxgene_interaction_models.Rmd",
    gx_annot="data/metadata/01_tx_annot.tsv",
    coxdata="data/Rds/12_coxdata.Rds",
    tools="src/enrichment-analysis-functions.R",
    pw="results/diffex/07_pairwise_comparisons_allgenes.xlsx",
    ovr="results/diffex/08_one_vs_rest_allgenes.xlsx"
  output:
    html=expand("reports/13_{rep}.html", rep=interaction_models),
    csv=expand("results/survival/13_{rep}.csv", rep=interaction_models)
  shell:
    """
    Rscript {input.script}
    """

rule app_setup:
  input:
    "data/metadata/01_tx_annot.tsv",
    "data/external/ensembl_universal_ids_v94.txt",
    "results/survival/12_cox_allgenes.xlsx",
    "results/survival/13_multi_interaction_os.csv",
    "results/survival/13_multi_interaction_drs.csv",
    expand("results/diffex/{result}.xlsx", result=diffex_results),
    "data/Rds/12_coxdata.Rds"
  output:
    "shinyApp/VisualizePPBCgene/data/app_gx_annot.Rds",
    "shinyApp/VisualizePPBCgene/data/12_cox_allgenes.xlsx",
    "shinyApp/VisualizePPBCgene/data/13_multi_interaction_os.csv",
    "shinyApp/VisualizePPBCgene/data/13_multi_interaction_drs.csv",
    "shinyApp/VisualizePPBCgene/data/app_diffex_res_list.Rds",
    "shinyApp/VisualizePPBCgene/data/app_survival_sample_data.Rds",
    "shinyApp/VisualizePPBCgene/data/app_ensembl_tmmnorm_genesxsample.Rds",
    "shinyApp/VisualizePPBCgene/data/app_symbol_tmmnorm_genesxsample.Rds"
  shell:
    """
    Rscript shinyApp/VisualizePPBCgene/data_setup.R
    """

subgroups=["Basal", "Her2", "LumA", "LumB"]
sub_diffex=[
  "ppbc_inv_vs_non_prbc",
  "ppbc_inv_vs_prbc",
  "ppbc_inv_vs_rest"
]

rds_dir="/DATA/share/postpartumbc/data/Rds"

rule subgroup_diffex:
  input:
    dds="data/Rds/08_dds_ovr_inv_vs_rest.Rds",
    script="src/14_subgroup_diffex.R"
  output:
    dds_res=expand("{dir}/subgroup_diffex/14_dds_{subgroup}_{comp}.Rds", subgroup=subgroups, comp=sub_diffex, dir=rds_dir),
    ape_res=expand("{dir}/subgroup_diffex/14_ape_{subgroup}_{comp}.Rds", subgroup=subgroups, comp=sub_diffex, dir=rds_dir)
  shell:
    """
    export OMP_NUM_THREADS=1 
    Rscript {input.script}
    """

diffex_dir="/DATA/share/postpartumbc/results/diffex"
    
rule report_subgroup_comp:
  input:
    dds_res=expand("{dir}/subgroup_diffex/14_dds_{subgroup}_{comp}.Rds", subgroup=subgroups, comp=sub_diffex, dir=rds_dir),
    ape_res=expand("{dir}/subgroup_diffex/14_ape_{subgroup}_{comp}.Rds", subgroup=subgroups, comp=sub_diffex, dir=rds_dir),
    rmd="reports/14_subgroup_diffex_by_comparison.Rmd",
    script="src/14_batch_subgroup_diffex_reports.R"
  output:
    allgenes=expand("{dir}/14_subgroup_diffex_{comp}_allgenes.xlsx", comp=sub_diffex, dir=diffex_dir),
    siggenes=expand("{dir}/14_subgroup_diffex_{comp}_sig_genes.xlsx", comp=sub_diffex, dir=diffex_dir),
    reports=expand("/DATA/share/postpartumbc/reports/14_subgroup_diffex_{comp}.html",comp=sub_diffex)
  params:
    comp=[
      "ppbc_inv_vs_non_prbc",
      "ppbc_inv_vs_prbc",
      "ppbc_inv_vs_rest"
      ]
  shell:
    """
    Rscript {input.script}
    """
#Should work but doesn't
#Rscript -e rmarkdown::render(input = "{input.rmd}", output_file = "/DATA/share/postpartumbc/reports/14_subgroup_diffex_{params.comp}.html", params = params.comp)


rule cox_glm:
  input:
    coxdata="data/Rds/12_coxdata.Rds",
    script="src/15_cox_elastic_net.R"
  output:
    "data/Rds/15_glm_os_1000.Rds",
    "data/Rds/15_glm_os_5000.Rds",
    "data/Rds/15_glm_drs_1000.Rds",
    "data/Rds/15_glm_drs_5000.Rds"
  shell:
    """
    Rscript {input.script}
    """

#As above, but include only involution samples
rule inv_cox_glm:
  input:
    coxdata="data/Rds/12_invdata.Rds",
    script="src/15b_inv_cox_elastic_net.R"
  output:
    "data/Rds/15b_inv_glm_os_1000.Rds",
    "data/Rds/15b_inv_glm_os_5000.Rds",
    "data/Rds/15b_inv_glm_drs_1000.Rds",
    "data/Rds/15b_inv_glm_drs_5000.Rds"
  shell:
    """
    Rscript {input.script}
    """

#As above, but use all genes as input
rule all_inv_cox_glm:
  input:
    coxdata="data/Rds/12_invdata.Rds",
    script="src/15c_inv_cox_elastic_net.R"
  output:
    "data/Rds/15c_inv_glm_os_all.Rds",
    "data/Rds/15c_inv_glm_drs_all.Rds"
  shell:
    """
    Rscript {input.script}
    """

rule report_cox_glm:
  input:
    "data/Rds/15_glm_os_1000.Rds",
    "data/Rds/15_glm_os_5000.Rds",
    "data/Rds/15_glm_drs_1000.Rds",
    "data/Rds/15_glm_drs_5000.Rds",
    diffex_results=expand("results/diffex/{result}.xlsx", result=diffex_results),
    gw_surv=expand("results/survival/12_{rep}.csv", rep=genewise_cox),
    int_surv=expand("results/survival/13_{rep}.csv", rep=interaction_models),
    gx_annot="data/metadata/01_tx_annot.tsv",
    rmd="reports/15_cox_elastic_net.Rmd",
    script="src/rmarkdown.R"
  output:
    "results/survival/15_elastic_cox_features.xlsx",
    html="reports/15_cox_elastic_net.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    
rule report_inv_cox_glm:
  input:
    "data/Rds/15b_inv_glm_os_1000.Rds",
    "data/Rds/15b_inv_glm_os_5000.Rds",
    "data/Rds/15b_inv_glm_drs_1000.Rds",
    "data/Rds/15b_inv_glm_drs_5000.Rds",
    diffex_results=expand("results/diffex/{result}.xlsx", result=diffex_results),
    gw_surv=expand("results/survival/12_{rep}.csv", rep=genewise_cox),
    int_surv=expand("results/survival/13_{rep}.csv", rep=interaction_models),
    gx_annot="data/metadata/01_tx_annot.tsv",
    rmd="reports/15b_inv_cox_elastic_net.Rmd",
    script="src/rmarkdown.R"
  output:
    "results/survival/15b_inv_elastic_cox_features.xlsx",
    html="reports/15b_inv_cox_elastic_net.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"    

rule report_all_inv_cox_glm:
  input:
    "data/Rds/15c_inv_glm_os_all.Rds",
    "data/Rds/15c_inv_glm_drs_all.Rds",
    diffex_results=expand("results/diffex/{result}.xlsx", result=diffex_results),
    gw_surv=expand("results/survival/12_{rep}.csv", rep=genewise_cox),
    int_surv=expand("results/survival/13_{rep}.csv", rep=interaction_models),
    gx_annot="data/metadata/01_tx_annot.tsv",
    rmd="reports/15c_inv_cox_elastic_net.Rmd",
    script="src/rmarkdown.R"
  output:
    "results/survival/15c_all_inv_elastic_cox_features.xlsx",
    html="reports/15c_all_inv_cox_elastic_net.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}" 
    
enet_features=[
  "15_elastic_cox_features",
  "15b_inv_elastic_cox_features",
  "15c_all_inv_elastic_cox_features"
]

inv_bf = [
  "breastfeeding",
  "involution"
]

rule diffex_involution:
  input:
    script="src/rmarkdown.R",
    rmd="reports/09_diffex_time_involution.Rmd",
    dds="data/Rds/08_dds_ovr_inv_vs_rest.Rds",
    vsd="data/Rds/08_vsd_ovr.Rds",
    gx_annot="data/metadata/01_tx_annot.tsv",
    sets=expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    cp="data/Rds/color_palettes.Rds",
    sp="data/Rds/survival_colors.Rds"
  output:
    html="reports/09_diffex_time_involution.html",
    dds="data/Rds/09_dds_involution_duration.Rds",
    ape="data/Rds/09_ape_involution_duration.Rds",
    genes="results/diffex/09_diffex_involution_duration.xlsx",
    heatmap="results/diffex/figs/09_involution_duration/09_hm_involution_duration.pdf",
    volcano="results/diffex/figs/09_involution_duration/09_volcano_involution_duration.jpeg"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    
rule diffex_breastfeeding:
  input:
    script="src/rmarkdown.R",
    rmd="reports/09_diffex_breastfeeding_duration.Rmd",
    dds="data/Rds/08_dds_ovr_inv_vs_rest.Rds",
    vsd="data/Rds/08_vsd_ovr.Rds",
    gx_annot="data/metadata/01_tx_annot.tsv",
    sets=expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    cp="data/Rds/color_palettes.Rds",
    sp="data/Rds/survival_colors.Rds"
  output:
    html="reports/09_diffex_time_breastfeeding.html",
    dds="data/Rds/09_dds_breastfeeding_duration.Rds",
    ape="data/Rds/09_ape_breastfeeding_duration.Rds",
    genes="results/diffex/09_diffex_breastfeeding_duration.xlsx",
    heatmap="results/diffex/figs/09_breastfeeding_duration/09_hm_breastfeeding_duration.pdf",
    volcano="results/diffex/figs/09_breastfeeding_duration/09_volcano_breastfeeding_duration.jpeg"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"    

rule gene_unity_setup:
  input:
    script="src/rmarkdown.R",
    rmd="reports/16_gene_unity_setup.Rmd",
    diffex_results=expand("results/diffex/{result}.xlsx", result=diffex_results),
    invbf_diffex=expand("results/diffex/09_diffex_{result}_duration.xlsx", result=inv_bf),
    #gw_surv=expand("results/survival/12_{rep}.csv", rep=genewise_cox),
    gw_surv="results/survival/12_cox_allgenes.xlsx",
    int_surv=expand("results/survival/13_{rep}.csv", rep=interaction_models),
    efeat=expand("results/survival/{feat}.xlsx", feat=enet_features),
    subdiffex=expand("{dir}/14_subgroup_diffex_{comp}_allgenes.xlsx", comp=sub_diffex, dir=diffex_dir),
    gx_annot="data/metadata/01_tx_annot.tsv",
    coxdata="data/Rds/12_coxdata.Rds",
    dds="data/Rds/08_dds_ovr_inv_vs_rest.Rds"
  output:
    html="reports/16_gene_unity_setup.html",
    aggdata="data/Rds/16_gene_report_environment.RData"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"   


rule gene_reports:
  input:
    aggdata="data/Rds/16_gene_report_environment.RData",
    script="src/17_batch_gene_reports.R",
    genes=ancient("reports/genes_to_report")
  #output:
  #  directory("reports/gene_reports")
  shell:
    "Rscript {input.script}"
    
rule trust_setup:
  output:
    bcrtcrfa="bin/TRUST4/hg38_bcrtcr.fa",
    imgt="bin/TRUST4/human_IMGT+C.fa",
    trust="bin/TRUST4/run-trust4"
  shell:
    """
    git clone https://github.com/liulab-dfci/TRUST4.git
    mkdir -p bin
    mv -v TRUST4 ./bin
    """

rule trust:
  input:
    fq=expand("data/RAW/{sample}.fastq.gz", sample=config["samples"]),
    bcrtcrfa=ancient("bin/TRUST4/hg38_bcrtcr.fa"),
    imgt=ancient("bin/TRUST4/human_IMGT+C.fa"),
    trust=ancient("bin/TRUST4/run-trust4")
  params:
    threads=16,
    outdir="data/TRUST/"
  shell:
    """
    mkdir -p {params.outdir} &&
     #for f in `ls $PROJDIR/data/RAW/*.R1.fastq.gz`
    for f in `ls data/RAW/*.R1.fastq.gz`; do basefile="$(basename -- $f)"; {input.trust} -u $f -t {params.threads} -f {input.bcrtcrfa} --ref {input.imgt} -o {params.outdir}$basefile; done
    """

rule trust_report:
  input:
    #directory("data/TRUST"),
    fqdata="data/metadata/01_sample_annot.tsv",
    sampledata="data/metadata/05_sample_annot_filtered.csv",
    survdata="data/Rds/04_survdata.Rds",
    preexcluded_samples="data/metadata/01_pre_excluded_samples.csv",
    discarded_samples="data/metadata/02_discarded_samples.csv",
    ihc_outliers="data/metadata/03_removed_pam50_outliers.csv",
    rmd="reports/18_BCR_clonality.Rmd",
    script="src/rmarkdown.R"
  output:
    alltrust="data/Rds/18_alltrust.Rds",
    trustdata="data/Rds/18_trustdata.Rds",
    trustexcel="results/TRUST/18_TRUST_results.xlsx",
    html="reports/18_BCR_clonality.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"  

#Does not work as a rule but does work in shell directly, fixme    
rule workflow_diagram:
  conda:
    "envs/environment.yml"
  shell:
    "snakemake --use-conda --dag | dot -Tsvg > dag.svg"
    #To view via shell: display dag.svg
    

### Spatial analyses ###

# The spatial analysis gets its own subdirectory in ./reports and within the data (sub)directories
# Raw data is in ./data/vectra

# Create a manageable data structure
rule organize_vectra:
  input:
    #directory("data/vectra/Batch 3 Leuven (Batch 1 in HALO)"), #Raw data export
    metadata="data/metadata/spatial/PPBC_sample_set.xlsx",
    rmd="reports/spatial/00_organize_vectra_samples.Rmd",
    script="src/rmarkdown.R"
  output:
    #dir("data/vectra/symlinks"),
    filedict="data/metadata/spatial/00_file_location_dictionary.csv",
    meta="data/metadata/spatial/00_vectra_metadata.csv",
    html="reports/spatial/00_organize_vectra_samples.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

vectra_table_pth = Path("data/metadata/spatial/00_file_location_dictionary.csv")
if vectra_table_pth.exists():
    vectra_table = pd.read_csv(vectra_table_pth)
    t_numbers = list(pd.unique(vectra_table['t_number']))
else:
    t_numbers = []

#QC checks based on aggregated summary files
rule summary_qc:
  input:
    filedict="data/metadata/spatial/00_file_location_dictionary.csv",
    rmd="reports/spatial/01_summary_QC.Rmd",
    script="src/rmarkdown.R"
  output:
    mpif26_summary="data/vectra/interim/summaries/01_MPIF26_batch1_summaries.csv",
    mpif27_summary="data/vectra/interim/summaries/01_MPIF27_batch1_summaries.csv",
    html="reports/spatial/01_summary_QC.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"

all_objects = expand(
  "data/vectra/interim/objects/{t_number}_{panel}_batch1.nc",
  t_number = t_numbers,
  panel = ['MPIF26', 'MPIF27'])

rule object_QC:
  input:
    rmd="reports/spatial/02_object_QC.Rmd",
    script="src/rmarkdown.R",
    mpif26_objects="data/vectra/interim/MPIF26_df.Rds",
    mpif27_objects="data/vectra/interim/MPIF27_df.Rds",
  output:
    html="reports/spatial/02_object_QC.html",
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"

rule all_objects:
  input: all_objects

rule load_objects:
  input:
    script = "src/spatial/import_halo_objects.py",
    objects="data/vectra/symlinks/objects/{t_number}_{panel}_{batch}_object_results.csv",
    summary="data/vectra/symlinks/summary/{t_number}_{panel}_{batch}_summary_results.csv",
  output:
    objects="data/vectra/interim/objects/{t_number}_{panel}_{batch}.nc",
  shell:
    "export OMP_NUM_THREADS=1\n"
    "python {input.script} {input.summary} {input.objects} {output.objects}"

rule cell_type_density:
  input:
    script="src/spatial/model_density.R",
    objects="data/vectra/processed/objects_{panel}.Rds",
  output:
    tsv="results/spatial/density/{panel}.tsv",
  shell:
    "mkdir -p results/spatial/density\n"
    "Rscript {input.script} {input.objects} {output}"

rule count_cells:
  input:
    all_objects,
    "src/spatial/read_cells.R",
    script="src/spatial/counts_cells.R",
  output:
    csv="results/spatial/cell_counts_by_marker.csv",
  shell:
    "mkdir -p results/spatial/\n"
    "Rscript {input.script} data/vectra/interim/objects/ {output.csv}"

rule convert_objects_nc_to_rds:
  input:
    object_nc=expand(
      "data/vectra/interim/objects/{t_number}_{{panel}}_batch1.nc",
      t_number = t_numbers),
    script="src/spatial/read_cells.R"
  output:
    "data/vectra/interim/{panel}_df.Rds",
  shell:
    "Rscript {input.script} {input.object_nc} {output}"

rule report_marker_correction:
  input:
    mpif26="data/vectra/interim/MPIF26_df.Rds",
    mpif27="data/vectra/interim/MPIF27_df.Rds",
    cell_count_by_marker="results/spatial/cell_counts_by_marker.csv",
    rmd="reports/spatial/03a_marker_correction.Rmd",
    script="src/rmarkdown.R"
  output:
    html="reports/spatial/03a_marker_correction.html",
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"

rule test_marker_correction_results:
  input:
    mpif26="data/vectra/processed/objects_MPIF26.Rds",
    mpif27="data/vectra/processed/objects_MPIF27.Rds",
    rmd="reports/spatial/03b_test_marker_correction.Rmd",
    script="src/rmarkdown.R"
  output:
    html="reports/spatial/03b_test_marker_correction.html",
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"

rule call_cell_types:
  input:
    script="src/spatial/call_cell_types.R",
    objects="data/vectra/interim/{panel}_df.Rds",
  output:
    "data/vectra/processed/objects_{panel}.Rds",
  shell:
    "Rscript {input.script} {input.objects} {wildcards.panel} {output}"
    
rule define_cell_types_report:
  input:
    mpif26="data/vectra/processed/objects_MPIF26.Rds",
    mpif27="data/vectra/processed/objects_MPIF27.Rds",
    rmd="reports/spatial/04_define_cell_types.Rmd",
    script="src/rmarkdown.R"
  output:
    html="reports/spatial/04_define_cell_types.html"
  shell:
    "Rscript {input.script} {input.rmd} $(realpath -s {output.html})"    
