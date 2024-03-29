#### Set up ####

# Simplify the fastq dir structure
rule symlinks:
  input:
    rmd="reports/rnaseq/00_symlink_fastqs.Rmd",
    script="src/utils/rmarkdown.R",
    sample_paths="data/rnaseq/metadata/sample_paths.txt"
  output:
    html="reports/rnaseq/00_symlink_fastqs.html"
  params:
    raw_dir="data/rnaseq/RAW",
    sym_dir="data/rnaseq/fastq"  
  shell:
   "mkdir -p data/rnaseq/fastq\n"
   "Rscript {input.script} {input.rmd} $PWD/{output.html}"
   " --sample_paths '{input.sample_paths}'"
   " --raw_dir '{params.raw_dir}'"
   " --sym_dir '{params.sym_dir}'"

# Config file for salmon quantification   
rule create_config:
  input:
    script="src/rnaseq/generate_config.py"
  output:
    config="config.yaml"
  shell:
    "python3 {input.script}"

#### Salmon quantification ####

rule salmon_index:
  input:
    script="src/rnaseq/build_salmon_index.sh"
  output:
    directory("data/external/index/grch38_index")
  conda:
    # File (or softlink) must be in the same dir as rnaseq.smk.py
    "salmon.yml"
  shell:
    "bash {input.script}"

rule salmon_quant:
  input:
    fastq = "data/rnaseq/fastq/{sample}.fastq.gz",
    index="data/external/index/grch38_index/"
  output:
    "data/rnaseq/salmon/{sample}/quant.sf"
  params:
    outdir = "data/rnaseq/salmon/{sample}"
  conda:
    "salmon.yml"
  threads: 16
  shell:
     "salmon quant -i {input.index} -l A --gcBias --seqBias "
     "--validateMappings --posBias -r {input.fastq} -o {params.outdir}"
     
#### Metadata and tximport ####
      
rule rna_metadata:
  input:
    sampledata="data/external/sample_data.tsv",
    patientdata="data/external/patient_data.tsv",
    rmd="reports/rnaseq/01_rna_metadata.Rmd",
    script="src/utils/rmarkdown.R",
    read_rna_sd="src/rnaseq/read_rna_sampledata.R",
    read_patient_data="src/utils/read_patient_data.R",
    pre_excluded_samples="data/external/pre_excluded_samples.csv"
  params:
    salmondir="data/rnaseq/salmon"
  output:
    html="reports/rnaseq/01_rna_metadata.html",
    rnaMeta="data/rnaseq/metadata/01_rnaMeta.Rds"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --sampledata '{input.sampledata}'"
    " --patientdata '{input.patientdata}'"
    " --salmondir '{params.salmondir}'"
    " --pre_excluded_samples '{input.pre_excluded_samples}'"

rule tximport:
  input:
    expand("data/rnaseq/salmon/{sample}/quant.sf", sample=config['samples']),
    rnaMeta="data/rnaseq/metadata/01_rnaMeta.Rds",
    biomart="data/external/gene_ref/biomart_ensemblgenes_v94.txt",
    rmd="reports/rnaseq/01b_tximport.Rmd",
    script="src/utils/rmarkdown.R"
  params:
    salmondir="data/rnaseq/salmon"
  output:
    gene_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    tx="data/rnaseq/interim/01_tx.Rds",
    html="reports/rnaseq/01b_tximport.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"    
    " --rnaMeta '{input.rnaMeta}'"
    " --salmondir '{params.salmondir}'"
    " --biomart '{input.biomart}'"
    
#### Quality control ####

rule fastqc:
  input:
    fastq="data/rnaseq/fastq/{sample}.fastq.gz"
  output:
    html = "results/rnaseq/fastqc/{sample}_fastqc.html"
  threads: 16
  conda:
    "fastqc.yml"         
  shell:
    """
    mkdir -p results/rnaseq/fastqc; fastqc {input.fastq} -o results/rnaseq/fastqc/ -t {threads}
    """

rule multiqc:
  input:
    fastqc=expand("results/rnaseq/fastqc/{sample}_fastqc.html", sample=config['samples'])
  params:
    fastqc_dir="results/rnaseq/fastqc/",
    salmon_dir = "data/rnaseq/salmon/"
  output:
    html="results/rnaseq/multiqc_report.html",
    alignstats="results/rnaseq/multiqc_report_data/multiqc_general_stats.txt",
    or_summary="results/rnaseq/multiqc_report_data/mqc_fastqc_overrepresented_sequencesi_plot_1.txt" 
  conda:
    "multiqc.yml"             
  shell:
    # -o sets the results directory, but be careful when using it with {output}
    # If passed {output.html}, Snakemake will search for results/rnaseq/multiqc.html,
    # but the actual file will be placed in results/rnaseq/results/rnaseq/multiqc.html
    """multiqc {params.fastqc_dir} {params.salmon_dir} -o results/rnaseq -n multiqc_report.html --force"""

rule fastqcr:
  input:
    script="src/rnaseq/02_fastqcr_aggregate.R",
  output:
    fastqcr="data/rnaseq/interim/02_fastqcr.Rds"
  shell:
    "Rscript {input.script}"
    
rule QC_salmon:
  input:
    script="src/utils/rmarkdown.R",
    rmd="reports/rnaseq/02_QC_salmon.Rmd",
    tx="data/rnaseq/interim/01_tx.Rds", 
    rnaMeta="data/rnaseq/metadata/01_rnaMeta.Rds",
    alignstats="results/rnaseq/multiqc_report_data/multiqc_general_stats.txt",
    or_sum="results/rnaseq/multiqc_report_data/mqc_fastqc_overrepresented_sequences_plot_1.txt",
    fastqcr="data/rnaseq/interim/02_fastqcr.Rds",
    lib="src/rnaseq/fastqcr_utils.R",
    pre_excluded_samples="data/external/pre_excluded_samples.csv"    
  output:
    discarded="data/rnaseq/metadata/02_discarded_samples.csv",
    meta_filtered = "data/rnaseq/interim/02_sample_annot_filtered.Rds",
    tx_clean= "data/rnaseq/interim/02_tx_clean.Rds",
    html="reports/rnaseq/02_QC_salmon.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --tx {input.tx}"
    " --rnaMeta {input.rnaMeta}"
    " --alignstats {input.alignstats}"
    " --or_sum '{input.or_sum}'"
    " --fastqcr '{input.fastqcr}'"

#### DeseqDataset and PAM50 ####
    
rule create_dds:
  input:
    rmd="reports/rnaseq/02b_create_dds.Rmd",
    tx_clean= "data/rnaseq/interim/02_tx_clean.Rds",
    meta_filtered = "data/rnaseq/interim/02_sample_annot_filtered.Rds",
    script="src/utils/rmarkdown.R"
  output:
    dds_qc = "data/rnaseq/interim/02_QC_dds.Rds",
    html="reports/rnaseq/02b_create_dds.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --tx {input.tx_clean}"
    " --meta {input.meta_filtered}"

rule pam50:
  input:
    script="src/utils/rmarkdown.R",
    rmd="reports/rnaseq/03_PAM50_IHC.Rmd",
    dds="data/rnaseq/interim/02_QC_dds.Rds",
    entrez="data/external/gene_ref/entrez_biomart_ensemblgenes_v94.txt"
  output:
    # Samples discarded due to clear incongruity
    removed="data/rnaseq/metadata/03_removed_pam50_outliers.csv",
    # New dds post Pam50 determination
    dds="data/rnaseq/interim/03_dds_PAM50.Rds",
    # New metadata post PAM50 determination
    pamsd="data/rnaseq/metadata/03_sample_annot_filtered_PAM50.csv",
    html="reports/rnaseq/03_PAM50_IHC.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --dds {input.dds}"
    " --entrez {input.entrez}"

#### Color configuration ####

rule color_palettes:
  input:
    rmd="reports/rnaseq/03b_color_palettes.Rmd",
    script="src/utils/rmarkdown.R",
    dds="data/rnaseq/interim/03_dds_PAM50.Rds"
  output:
    html="reports/rnaseq/03b_color_palettes.html",
    color_palettes="data/rnaseq/interim/color_palettes.Rds"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --dds {input.dds}"

rule survival_colors:
  input:
    script="src/rnaseq/survival_colors.R"
  output:
    colors="data/rnaseq/interim/survival_colors.Rds"
  shell:
    "Rscript {input.script}"
    
#### Survival and ESTIMATE ####

rule estimate:
  input:
    script="src/utils/rmarkdown.R",
    rmd="reports/rnaseq/04a_estimate_tumor_purity.Rmd",
    dds="data/rnaseq/interim/03_dds_PAM50.Rds",
    genes="data/rnaseq/metadata/01_gene_annot.tsv",
    color_palette="data/rnaseq/interim/color_palettes.Rds",
  output:
    hugo="data/rnaseq/interim/hugo_fpkm.txt",
    estimate_filter="results/rnaseq/ESTIMATE/filterCommonGenes.gct",
    estimate_score="results/rnaseq/ESTIMATE/results_estimate_score.gct",
    processed_scores="data/rnaseq/interim/04a_estimate_scores.Rds",
    html="reports/rnaseq/04a_estimate_tumor_purity.html",
    dds="data/rnaseq/interim/04_dds_PAM50_est.Rds"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --dds {input.dds}"
    " --genes {input.genes}"
    " --color_palette {input.color_palette}"

rule uni_survival:
  input:
    script="src/utils/rmarkdown.R",
    rmd="reports/rnaseq/04b_univariate_survival.Rmd",
    dds="data/rnaseq/interim/04_dds_PAM50_est.Rds",
    scores="data/rnaseq/interim/04a_estimate_scores.Rds",
    color_palette="data/rnaseq/interim/color_palettes.Rds",
    survival_colors="data/rnaseq/interim/survival_colors.Rds"
  output:
    excluded="data/rnaseq/metadata/04_samples_excluded_survival.csv",
    coxdata="data/rnaseq/interim/04_survdata.Rds",
    html="reports/rnaseq/04b_univariate_survival.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}" 
    " --dds {input.dds}"
    " --scores {input.scores}"
    " --color_palette {input.color_palette}"
    " --survival_colors {input.survival_colors}"

rule multi_survival:
  input:
    script="src/utils/rmarkdown.R",
    rmd="reports/rnaseq/04c_multivariate_survival.Rmd",
    coxdata="data/rnaseq/interim/04_survdata.Rds",
    color_palette="data/rnaseq/interim/color_palettes.Rds",
    dds="data/rnaseq/interim/04_dds_PAM50_est.Rds",
    survival_colors="data/rnaseq/interim/survival_colors.Rds"
  output:
    html="reports/rnaseq/04c_multivariate_survival.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --coxdata {input.coxdata}"
    " --dds {input.dds}"
    " --color_palette {input.color_palette}"
    " --survival_colors {input.survival_colors}"

#### Batch effects ####
    
rule batch_effects:
  input:
    script="src/utils/rmarkdown.R",
    rmd="reports/rnaseq/05_batch_effects.Rmd",
    dds="data/rnaseq/interim/04_dds_PAM50_est.Rds",
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    sampledata="data/external/sample_data.tsv",
    cp="data/rnaseq/interim/color_palettes.Rds"
  output:
    metadata="data/rnaseq/metadata/05_sample_annot_filtered.csv",
    dds="data/rnaseq/interim/05_dds_PAM50_batch.Rds",
    html="reports/rnaseq/05_batch_effects.html",
    pcas="results/rnaseq/dimensionality/pca_rnaseq.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --dds {input.dds}"
    " --gx_annot {input.gx_annot}"
    " --cp {input.cp}"
    " --sampledata {input.sampledata}"
    " --pcas {output.pcas}"

#### Differential expression ####

rule min_count_threshold:
  input:
    script = "src/utils/rmarkdown.R",
    rmd = "reports/rnaseq/05b_minimum_count_threshold.Rmd",
    dds = "data/rnaseq/interim/05_dds_PAM50_batch.Rds"
  output:
    dds = "data/rnaseq/interim/05b_dds_filtered.Rds",
    html = "reports/rnaseq/05b_minimum_count_threshold.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
     " --dds {input.dds}"

# Likelihood ratio tests
# Identify DEGs that are differentially expressed in at least one group

# Setting the number of threads with OMP_NUM_THREADS prevents unwanted parallelization    
rule lrt_diffex:
  input:
    dds="data/rnaseq/interim/05b_dds_filtered.Rds",
    gene_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    script = "src/rnaseq/diffex_lrt.R"
  output:
    "data/rnaseq/interim/06_vsd.Rds",
    "data/rnaseq/interim/06_vsd_nolac.Rds",
    "data/rnaseq/interim/06_ddsLRT.Rds",
    "data/rnaseq/interim/06_ddsLRT_nolac.Rds",
    "data/rnaseq/interim/06_ddsLRT_pam.Rds",
    "data/rnaseq/interim/06_ddsLRT_pam_nolac.Rds",
    "data/rnaseq/interim/06_ddsLRT_batch.Rds",
    "data/rnaseq/interim/06_ddsLRT_batch_nolac.Rds"
  shell:
    """
    export OMP_NUM_THREADS=1 
    Rscript {input.script}
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

# Report the LRT results
# This notebook has serious bloat but no time to refactor it now
rule lrt_report:
  input:
    expand("data/rnaseq/interim/06_{dds}.Rds", dds=dds_lrt),
    expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    dds="data/rnaseq/interim/05b_dds_filtered.Rds",
    vsd="data/rnaseq/interim/06_vsd.Rds",
    #ImmPort database: https://www.innatedb.com/redirect.do?go=resourcesGeneLists
    immune_genes="data/external/gene_ref/InnateDB_genes.csv",
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    cp="data/rnaseq/interim/color_palettes.Rds",
    sp="data/rnaseq/interim/survival_colors.Rds",
    tools="src/rnaseq/deseq_report_functions.R",
    rmd="reports/rnaseq/06_diffex_lrt.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    sig_genes="results/rnaseq/diffex/06_LRT_sig_genes.xlsx",
    all_genes="results/rnaseq/diffex/06_LRT_allgenes.xlsx",
    pathways="results/rnaseq/diffex/06_LRT_pathways.xlsx",
    html="reports/rnaseq/06_diffex_lrt.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --dds {input.dds}"
    " --gx_annot {input.gx_annot}"
    " --immune_genes {input.immune_genes}"
    " --cp {input.cp}"
    " --sp {input.sp}"
    " --tools {input.tools}"
    " --vsd {input.vsd}"


# Pairwise comparisons between all groups
pairwise_refs = ["npbc", "ppbcdl", "prbc"]

pairwise_apeglm = [
  "prbc_vs_npbc", "ppbcdl_vs_npbc", "ppbcpw_vs_npbc", 
  "prbc_vs_ppbcdl", "ppbcpw_vs_prbc", "ppbcpw_vs_ppbcdl"
]
    
rule diffex_pairwise:
  input:
    dds="data/rnaseq/interim/05b_dds_filtered.Rds",
    script="src/rnaseq/diffex_pairwise.R"
  output:
    pairwise_dds=expand("data/rnaseq/interim/07_dds_pairwise_ref_{pw}.Rds", pw=pairwise_refs),
    apeglm_results=expand("data/rnaseq/interim/07_ape_{pw}.Rds", pw=pairwise_apeglm)
  shell:
    """
    export OMP_NUM_THREADS=1 
    Rscript {input.script}
    """
    
pairwise_prefixes = [
  "inv_lac", "inv_nonprbc", "inv_prbc",
  "lac_nonprbc", "lac_prbc", "prbc_nonprbc"
]

# Report on pairwise Wald tests
rule pairwise_report:
  input:
    pairwise_dds=expand("data/rnaseq/interim/07_dds_pairwise_ref_{pw}.Rds", pw=pairwise_refs),
    sets=expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    immune_genes="data/external/gene_ref/InnateDB_genes.csv",
    cp="data/rnaseq/interim/color_palettes.Rds",
    sp="data/rnaseq/interim/survival_colors.Rds",
    vsd="data/rnaseq/interim/06_vsd.Rds",
    tools="src/rnaseq/deseq_report_functions.R",
    apeglm_results=expand("data/rnaseq/interim/07_ape_{pw}.Rds", pw=pairwise_apeglm),
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    rmd="reports/rnaseq/07_diffex_pairwise.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    all_genes="results/rnaseq/diffex/07_pairwise_comparisons_allgenes.xlsx",
    sig_genes="results/rnaseq/diffex/07_pairwise_comparisons_sig_genes.xlsx",
    volcano=expand("results/rnaseq/diffex/volc_{pw}.pdf", pw=pairwise_prefixes),
    pathways="results/rnaseq/diffex/07_pairwise_comparisons_pathways.xlsx",
    html="reports/rnaseq/07_diffex_pairwise.html"
  shell:
     "Rscript {input.script} {input.rmd} $PWD/{output.html}"
     " --gx_annot {input.gx_annot}"
     " --immune_genes {input.immune_genes}"
     " --cp {input.cp}"
     " --sp {input.sp}"
     " --tools {input.tools}"
     " --vsd {input.vsd}"

# One vs rest comparisons
# a single study group is compared against all other samples pooled together

ovr_comps = ["inv_vs_rest", "prbc_vs_rest", "lac_vs_rest", "nonprbc_vs_rest"]
     
rule diffex_one_vs_rest:
  input:
    dds="data/rnaseq/interim/05b_dds_filtered.Rds",
    script="src/rnaseq/diffex_onevsrest.R"
  output:
    dds_ovr=expand("data/rnaseq/processed/08_dds_ovr_{comp}.Rds", comp=ovr_comps),
    ape_ovr=expand("data/rnaseq/processed/08_ape_ovr_{comp}.Rds", comp=ovr_comps)
  shell:
    """
    export OMP_NUM_THREADS=1 
    Rscript {input.script}
    """

# Report on Wald tests from one vs rest comparisons
    
ovr_reports = ["inv_rest", "prbc_rest", "lac_rest", "nonprbc_rest"]

rule one_vs_rest_report:
  input:
    ape_ovr=expand("data/rnaseq/processed/08_ape_ovr_{comp}.Rds", comp=ovr_comps),
    dds_ovr=expand("data/rnaseq/processed/08_dds_ovr_{comp}.Rds", comp=ovr_comps),
    sets=expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    immune_genes="data/external/gene_ref/InnateDB_genes.csv",
    cp="data/rnaseq/interim/color_palettes.Rds",
    sp="data/rnaseq/interim/survival_colors.Rds",
    vsd="data/rnaseq/interim/06_vsd.Rds",
    tools="src/rnaseq/deseq_report_functions.R",
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    rmd="reports/rnaseq/08_diffex_onevsrest.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    vsd_ovr="data/rnaseq/interim/08_vsd_ovr.Rds",
    allgenes="results/rnaseq/diffex/08_one_vs_rest_allgenes.xlsx",
    sig_genes="results/rnaseq/diffex/08_one_vs_rest_sig_genes.xlsx",
    pathways="results/rnaseq/diffex/08_one_vs_rest_pathways.xlsx",
    heatmaps=expand("results/rnaseq/diffex/hm_{ovr}.pdf", ovr=ovr_reports),
    volcano=expand("results/rnaseq/diffex/volc_{ovr}.pdf", ovr=ovr_reports),
    html="reports/rnaseq/08_diffex_onevsrest.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
     " --gx_annot {input.gx_annot}"
     " --immune_genes {input.immune_genes}"
     " --cp {input.cp}"
     " --sp {input.sp}"
     " --tools {input.tools}"
     " --vsd {input.vsd}"

# Examine whether the IG gene signature is correlated with lactation
rule ig_milk:
  input:
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    tools="src/rnaseq/deseq_report_functions.R",
    lrt="results/rnaseq/diffex/06_LRT_allgenes.xlsx",
    vsd="data/rnaseq/interim/08_vsd_ovr.Rds",
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    cp="data/rnaseq/interim/color_palettes.Rds",
    sp="data/rnaseq/interim/survival_colors.Rds",
    ivr="results/rnaseq/diffex/08_one_vs_rest_allgenes.xlsx",
    rmd="reports/rnaseq/08b_ig_milk_genes.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    milk_heatmap="results/rnaseq/diffex/milk_vs_IG_genes.pdf",
    normat="data/rnaseq/processed/08_sfnorm_genesymbol_countmatrix.Rds",
    ig_milk_genes="data/rnaseq/processed/08_ig_and_milk_genes.Rds",
    cor_ig_milk_all="results/rnaseq/diffex/cor_ig_milk_all.pdf",
    cor_ig_milk_inv="results/rnaseq/diffex/cor_ig_milk_inv.pdf",
    html="reports/rnaseq/08b_ig_milk_genes.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --dds {input.dds}"
    " --tools {input.tools}"
    " --lrt {input.lrt}"
    " --vsd {input.vsd}"
    " --gx_annot {input.gx_annot}"
    " --sp {input.sp}"    
    " --cp {input.cp}"
    " --ivr {input.ivr}"

subgroups=["Basal", "Her2", "LumA", "LumB"]
sub_diffex=[
  "ppbcpw_vs_npbc",
  "ppbcpw_vs_prbc",
  "ppbcpw_vs_rest"
]

# Diffex between PAM50 subgroups
# Examine whether (for example) basal samples are different between study groups
rule subgroup_diffex:
  input:
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    script="src/rnaseq/subgroup_diffex.R"
  output:
    dds_res=expand("{dir}/subgroup_diffex/dds_{subgroup}_{comp}.Rds",
      subgroup=subgroups, comp=sub_diffex, dir="data/rnaseq/processed/subgroup_diffex"),
    ape_res=expand("{dir}/subgroup_diffex/ape_{subgroup}_{comp}.Rds",
      subgroup=subgroups, comp=sub_diffex, dir="data/rnaseq/processed/subgroup_diffex")
  shell:
    """
    export OMP_NUM_THREADS=1 
    Rscript {input.script}
    """

# Report per sub_diffex.    
rule report_subgroup_comp:
  input:
    dds_res=expand("data/rnaseq/processed/subgroup_diffex/dds_{subgroup}_{comp}.Rds",
      subgroup=subgroups, comp=sub_diffex),
    ape_res=expand("data/rnaseq/processed/subgroup_diffex/ape_{subgroup}_{comp}.Rds",
      subgroup=subgroups, comp=sub_diffex),
    sets=expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    immune_genes="data/external/gene_ref/InnateDB_genes.csv",
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    cp="data/rnaseq/interim/color_palettes.Rds",
    sp="data/rnaseq/interim/survival_colors.Rds",
    lib="src/rnaseq/enrichment-analysis-functions.R",
    rmd="reports/rnaseq/14_subgroup_diffex_by_comparison.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    allgenes="results/rnaseq/diffex/14_subgroup_diffex_{comp}_allgenes.xlsx",
    siggenes="results/rnaseq/diffex/14_subgroup_diffex_{comp}_sig_genes.xlsx",
    html="reports/rnaseq/14_subgroup_diffex_{comp}.html",
    # Snakemake complains if not all output files have the same wildcard
    vp_basal="results/rnaseq/diffex/14_volc_{comp}_subgroup_Basal.pdf",
    vp_her2="results/rnaseq/diffex/14_volc_{comp}_subgroup_Her2.pdf",
    vp_lumA="results/rnaseq/diffex/14_volc_{comp}_subgroup_LumA.pdf",
    vp_lumB="results/rnaseq/diffex/14_volc_{comp}_subgroup_LumB.pdf"
  shell:
     "Rscript {input.script} {input.rmd} $PWD/{output.html}"
     " --comparison {wildcards.comp}"
     " --gx_annot {input.gx_annot}"
     " --immune_genes {input.immune_genes}"
     " --cp {input.cp}"
     " --sp {input.sp}"
     " --lib {input.lib}"

# Examine DEGs responsive to involution duration in PPBCpw
rule diffex_involution:
  input:
    script="src/utils/rmarkdown.R",
    rmd="reports/rnaseq/09_diffex_time_involution.Rmd",
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    vsd="data/rnaseq/interim/08_vsd_ovr.Rds",
    sets=expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    immune_genes="data/external/gene_ref/InnateDB_genes.csv",
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    cp="data/rnaseq/interim/color_palettes.Rds",
    sp="data/rnaseq/interim/survival_colors.Rds",
    lib="src/rnaseq/enrichment-analysis-functions.R"
  output:
    html="reports/rnaseq/09_diffex_time_involution.html",
    dds="data/rnaseq/processed/09_dds_involution_duration.Rds",
    ape="data/rnaseq/processed/09_ape_involution_duration.Rds",
    genes="results/rnaseq/diffex/09_diffex_involution_duration.xlsx",
    heatmap="results/rnaseq/diffex/09_hm_involution_duration.pdf",
    volcano="results/rnaseq/diffex/09_volcano_involution_duration.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --gx_annot {input.gx_annot}"
    " --immune_genes {input.immune_genes}"
    " --cp {input.cp}"
    " --sp {input.sp}"
    " --dds {input.dds}"
    " --vsd {input.vsd}"
    " --lib {input.lib}"

# Examine DEGs responsive to breastfeeding duration in PPBCpw
rule diffex_breastfeeding:
  input:
    script="src/utils/rmarkdown.R",
    rmd="reports/rnaseq/09b_diffex_breastfeeding_duration.Rmd",
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    vsd="data/rnaseq/interim/08_vsd_ovr.Rds",
    sets=expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    immune_genes="data/external/gene_ref/InnateDB_genes.csv",
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    cp="data/rnaseq/interim/color_palettes.Rds",
    sp="data/rnaseq/interim/survival_colors.Rds",
    lib="src/rnaseq/enrichment-analysis-functions.R"
  output:
    html="reports/rnaseq/09b_diffex_time_breastfeeding.html",
    dds="data/rnaseq/processed/09b_dds_breastfeeding_duration.Rds",
    ape="data/rnaseq/processed/09b_ape_breastfeeding_duration.Rds",
    genes="results/rnaseq/diffex/09b_diffex_breastfeeding_duration.xlsx",
    heatmap="results/rnaseq/diffex/09b_hm_breastfeeding_duration.pdf",
    volcano="results/rnaseq/diffex/09b_volcano_breastfeeding_duration.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --gx_annot {input.gx_annot}"
    " --immune_genes {input.immune_genes}"
    " --cp {input.cp}"
    " --sp {input.sp}"
    " --dds {input.dds}"
    " --vsd {input.vsd}"
    " --lib {input.lib}"

#### Pathway analysis ####

rule flexgsea:
  input:
    script = "src/rnaseq/deseq_flexgsea.R",
    lib = "src/rnaseq/deseq_report_functions.R",
    dds = "data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    gx_annot = "data/rnaseq/metadata/01_gene_annot.tsv",
    gene_sets = expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets)
  output:
    dir("results/rnaseq/flexgsea/deseq"),
    # Example results file
    "results/rnaseq/flexgsea/deseq/results/inv_vs_rest_canonpath_c2_results_flexdeseq.Rds",
    "results/rnaseq/flexgsea/deseq/input/geneSymbol_countmatrix.Rds"
  shell:
    """
    export OMP_NUM_THREADS=1 
    Rscript {input.script}
    """

rule flexgsea_report:
  input:
    # Takes a long time to rerun flexgsea with deseq and multiple comparisons
    ancient(dir("results/rnaseq/flexgsea/deseq")),
    ancient("results/rnaseq/flexgsea/deseq/results/inv_vs_rest_canonpath_c2_results_flexdeseq.Rds"),
    countMatrix = ancient("results/rnaseq/flexgsea/deseq/input/geneSymbol_countmatrix.Rds"),
    script = "src/utils/rmarkdown.R",
    rmd = "reports/rnaseq/09c_flexgsea_results.Rmd",
    gx_annot = "data/rnaseq/metadata/01_gene_annot.tsv"
  params:
    fdr_thresh = 0.25
  output:
    html = "reports/rnaseq/09c_flexgsea_results.html",
    overview = "results/rnaseq/flexgsea/deseq/flexgsea_aggregate_results.xlsx"
  shell:
     "Rscript {input.script} {input.rmd} $PWD/{output.html}"
     " --fdr_thresh '{params.fdr_thresh}'"
     " --gx_annot '{input.gx_annot}'"
     " --countMatrix '{input.countMatrix}'"

# Some exploration of TF enrichment
rule tfea:
  input:
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    one_vs_rest="results/rnaseq/diffex/08_one_vs_rest_allgenes.xlsx",
    sig_ovr="results/rnaseq/diffex/08_one_vs_rest_sig_genes.xlsx",
    rmd="reports/rnaseq/09d_cetf.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    #out="data/rnaseq/processed/cetf_out.Rds", #do not auto-delete on rerun
    html="reports/rnaseq/09d_cetf.html"
  shell:
    "export OMP_NUM_THREADS=1;"
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --dds '{input.dds}'"
    " --one_vs_rest '{input.one_vs_rest}'"
    " --sig_ovr '{input.sig_ovr}'"

rule gsea_tf:
  input:
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    ovr="results/rnaseq/diffex/08_one_vs_rest_allgenes.xlsx",
    gmt="data/external/gmt/c3.tft.v2022.1.Hs.symbols.gmt",
    gx_annot="data/rnaseq/processed/bx_annot.Rds",
    rmd="reports/rnaseq/09d_flexgsea_TFs.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    #outFile="data/rnaseq/processed/TF_flexgsea_invVsRest.Rds", #do not auto-delete on rerun
    html="reports/rnaseq/09d_flexgsea_TFs.html"
  shell:
    "export OMP_NUM_THREADS=1;"
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --dds '{input.dds}'"
    " --ovr '{input.ovr}'"
    " --gmt '{input.gmt}'"
    " --gx_annot '{input.gx_annot}'"

#### Cibersort Deconvolution ####

# Read in Cibersort results for plotting and survival analyses
# Large notebook: consider splitting
rule cibersortX:
  input:
    hugo="data/rnaseq/interim/hugo_fpkm.txt",
    #Results obtained by uploading hugo_fpkm.txt to the cibersort website with the enumerated parameters
    ciber_input="results/rnaseq/cibersortX/CIBERSORTx_Job6_Results.csv",
    abs_ciber_input="results/rnaseq/cibersortX/CIBERSORTx_Job7_Results.csv",
    coxdata="data/rnaseq/interim/04_survdata.Rds",
    metadata="data/rnaseq/metadata/05_sample_annot_filtered.csv",
    cp="data/rnaseq/interim/color_palettes.Rds",
    sp="data/rnaseq/interim/survival_colors.Rds",
    rmd="reports/rnaseq/10_CibersortX.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    report="reports/rnaseq/10_CibersortX.html",
    ciberdf="data/rnaseq/processed/10_ciberdf.Rds",
    abs_ciberdf="data/rnaseq/processed/10_abs_ciberdf.Rds",
    rel_hm="figures/supfigs/SupFig12a_cibersort_relative_heatmap.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --ciber_input '{input.ciber_input}'"
    " --abs_ciber_input '{input.abs_ciber_input}'"
    " --hugo '{input.hugo}'"
    " --coxdata '{input.coxdata}'"
    " --metadata '{input.metadata}'"
    " --cp '{input.cp}'"
    " --sp '{input.sp}'"
    " --rel_hm '{output.rel_hm}'"

#### Clustering of PPBCpw DEGs ####

#Only those results which include an inv vs something comparison
inv_comps = [
  "07_pairwise_comparisons_allgenes",
  "07_pairwise_comparisons_sig_genes",
  "08_one_vs_rest_allgenes",
  "08_one_vs_rest_sig_genes"
]

#Explore k-means clustering of differentially expressed genes and study group
#This notebook is also in need of splitting into smaller sections
rule inv_clustering:
  input:
    expand("results/rnaseq/diffex/{inv_comp}.xlsx", inv_comp=inv_comps),
    rmd="reports/rnaseq/11_clustering_involution.Rmd",
    script="src/utils/rmarkdown.R",
    cp="data/rnaseq/interim/color_palettes.Rds",
    sp="data/rnaseq/interim/survival_colors.Rds",
    tools="src/rnaseq/deseq_report_functions.R",
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    vsd="data/rnaseq/interim/08_vsd_ovr.Rds",
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    coxdata="data/rnaseq/interim/04_survdata.Rds",
    risktable="src/figures/faceted_risktable.R"
  output:
    cluster_heatmap="results/rnaseq/clustering/11_hm_clust_DEG_inv_vs_rest.pdf",
    DEGhm1="figures/Fig1a_DEG_heatmap.pdf",
    DEGhm2="figures/Fig1b_DEG_heatmap_kmeans.pdf",
    IGbarplot="figures/Fig1e_IG_cluster_barplot.pdf",
    IG_kaplan_PPBC_OS="figures/Fig1f_IG_cluster_PPBC_KM_OS.pdf",
    cluster_barplot="results/rnaseq/clustering/11_barplots_ig_clusters.pdf",
    cluster_results="results/rnaseq/clustering/11_inv_clusters.xlsx",
    cluster_survival="data/rnaseq/interim/11_ig_survdata.Rds",
    IG_kaplan_PPBC_DRS="figures/supfigs/Supfig9_IG_cluster_PPBC_KM_DRS.pdf",
    IG_kaplan_PAM50_OS="figures/supfigs/Supfig10a_IG_cluster_PAM50_KM_OS.pdf",
    IG_kaplan_PAM50_DRS="figures/supfigs/Supfig10b_IG_cluster_PAM50_KM_DRS.pdf",
    html="reports/rnaseq/11_clustering_involution.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --gx_annot {input.gx_annot}"
    " --cp {input.cp}"
    " --sp {input.sp}"
    " --dds {input.dds}"
    " --vsd {input.vsd}"
    " --tools {input.tools}"
    " --coxdata '{input.coxdata}'"
    " --risktable '{input.risktable}'"
    " --DEGhm1 '{output.DEGhm1}'"
    " --DEGhm2 '{output.DEGhm2}'"
    " --IGbarplot '{output.IGbarplot}'"
    " --IG_kaplan_PPBC_OS '{output.IG_kaplan_PPBC_OS}'"
    " --IG_kaplan_PPBC_DRS '{output.IG_kaplan_PPBC_DRS}'"
    " --IG_kaplan_PAM50_OS '{output.IG_kaplan_PAM50_OS}'"
    " --IG_kaplan_PAM50_DRS '{output.IG_kaplan_PAM50_DRS}'"
    
#### Survival analyses based on gene expression ####

# Add gene expression as columns to clinical covariates  
rule prepare_coxdata:
  input:
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    coxdata="data/rnaseq/interim/04_survdata.Rds",
    dds = "data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    script = "src/rnaseq/prepare_coxdata.R"
  output:
    coxdata="data/rnaseq/processed/12_coxdata.Rds", #includes gene expression
    invcoxdata="data/rnaseq/processed/12_invdata.Rds"
  shell:
    """
    Rscript {input.script}
    """

# Run the univariate and multivariate cox regressions
rule genewise_survival:
  input:
    coxdata="data/rnaseq/processed/12_coxdata.Rds",
    script = "src/rnaseq/genewise_survival.R"
  output:
    expand("data/rnaseq/processed/12_{res}.Rds",
      res=["multi_genewise_os", "uni_genewise_os", "multi_genewise_drs", "uni_genewise_drs"])
  shell:
    """
    Rscript {input.script}
    """

# As above, but with inv samples only
# Could be combined with genewise_survival with wildcards for brevity
rule inv_genewise_survival:
  input:
    coxdata="data/rnaseq/processed/12_invdata.Rds",
    lib = "src/rnaseq/genewise_survival.R",
    script = "src/rnaseq/inv_genewise_survival.R"
  output:
    expand("data/rnaseq/processed/12_{res}.Rds",
      res=["inv_multi_genewise_os", "inv_uni_genewise_os", "inv_multi_genewise_drs", "inv_uni_genewise_drs"])
  shell:
    """
    Rscript {input.script}
    """

genewise_cox =[
  "multi_genewise_os", "uni_genewise_os",
  "multi_genewise_drs", "uni_genewise_drs",
  "inv_multi_genewise_os", "inv_uni_genewise_os",
  "inv_multi_genewise_drs", "inv_uni_genewise_drs"
  ]

interaction_cox = [
  "uni_interaction_os",
  "uni_interaction_drs",
  "multi_interaction_os",
  "multi_interaction_drs"
]

# Interaction survival models, using the general formula
# outcome ~ gene*involution (univariate)
# outcome ~ clinical_covariates + gene*involution (multivariate)
# could be combined with the previous rules?
rule surv_inv_int:
  input:
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    coxdata="data/rnaseq/processed/12_coxdata.Rds",
    script="src/rnaseq/survival_involution_interaction.R"
  output:
    expand("data/rnaseq/processed/12_{m}.Rds", m=interaction_cox)
  shell:
    """
    Rscript {input.script}
    """

# Separately report the results of each genewise and interaction Cox models
rule report_genewise_survival:
  input:
    survival_results="data/rnaseq/processed/12_{cox}.Rds",
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    cp="data/rnaseq/interim/color_palettes.Rds",
    sp="data/rnaseq/interim/survival_colors.Rds",
    sets=expand("data/external/gmt/{gene_set}.gmt", gene_set=gene_sets),
    script="src/utils/rmarkdown.R",
    rmd="reports/rnaseq/12_genewise_survival.Rmd",
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    #either coxdata or invdata is loaded conditionally
    coxdata="data/rnaseq/processed/12_coxdata.Rds",
    invcoxdata="data/rnaseq/processed/12_invdata.Rds",
    tools="src/rnaseq/enrichment-analysis-functions.R",
    duplicates="src/rnaseq/summarize_duplicate_ids.R",
    pw="results/rnaseq/diffex/07_pairwise_comparisons_allgenes.xlsx",
    ovr="results/rnaseq/diffex/08_one_vs_rest_allgenes.xlsx"
  output:
    html="reports/rnaseq/12_{cox}.html",
    csv="results/rnaseq/survival/12_{cox}.csv",
    heatmap="results/rnaseq/survival/12_{cox}_heatmap.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --survival_results {input.survival_results}"
    " --gx_annot {input.gx_annot}"
    " --tools {input.tools}"
    " --duplicates {input.duplicates}"
    " --pw {input.pw}"
    " --ovr {input.ovr}"

# Combine all survival results and compare
rule aggregate_genewise_survival:
  input:
    csv=expand("results/rnaseq/survival/12_{rep}.csv", rep=genewise_cox+interaction_cox),
    script="src/utils/rmarkdown.R",
    rmd="reports/rnaseq/13_aggregate_genewise_survival.Rmd"
  params:
    resdir="results/rnaseq/survival"
  output:
    html="reports/rnaseq/13_aggregate_genewise_survival.html",
    coxres="results/rnaseq/survival/12_cox_allgenes.xlsx"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --resdir {params.resdir}"

#### Gene summary reports ####

# Add more identifier types for reporting
rule universal_ids:
  input:
    gx_annot="data/rnaseq/metadata/01_gene_annot.tsv",
    universal_ids="data/external/gene_ref/ensembl_universal_ids_v94.txt",
    script="src/rnaseq/universal_gene_ids.R"
  output:
    bx_annot="data/rnaseq/processed/bx_annot.Rds"
  shell:
    "Rscript {input.script}"

# Create an environment for generating gene reports
rule gene_unity_setup:
  input:
    script="src/utils/rmarkdown.R",
    rmd="reports/rnaseq/16_gene_unity_setup.Rmd",
    diffex_lrt="results/rnaseq/diffex/06_LRT_allgenes.xlsx",
    diffex_pw="results/rnaseq/diffex/07_pairwise_comparisons_allgenes.xlsx",
    diffex_ovr="results/rnaseq/diffex/08_one_vs_rest_allgenes.xlsx",
    inv_diffex="results/rnaseq/diffex/09_diffex_involution_duration.xlsx",
    bf_diffex="results/rnaseq/diffex/09b_diffex_breastfeeding_duration.xlsx",
    subdiffex=expand("results/rnaseq/diffex/14_subgroup_diffex_{comp}_allgenes.xlsx", 
      comp=sub_diffex),
    coxres="results/rnaseq/survival/12_cox_allgenes.xlsx",
    bx_annot="data/rnaseq/processed/bx_annot.Rds",
    coxdata="data/rnaseq/processed/12_coxdata.Rds",
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds"
  output:
    html="reports/rnaseq/16_gene_unity_setup.html",
    aggdata="data/rnaseq/processed/16_gene_report_environment.RData"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --dds {input.dds}"
    " --coxdata {input.coxdata}"
    " --gx_annot {input.bx_annot}"
    " --coxres {input.coxres}"
    " --diffex_lrt {input.diffex_lrt}"
    " --diffex_pw {input.diffex_pw}"
    " --diffex_ovr {input.diffex_ovr}"
    " --inv_diffex {input.inv_diffex}"
    " --bf_diffex {input.bf_diffex}"

rule gene_reports:
  input:
    aggdata="data/rnaseq/processed/16_gene_report_environment.RData",
    bx_annot="data/rnaseq/processed/bx_annot.Rds",
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    rmd="reports/rnaseq/17_gene_report_template.Rmd",
    script="src/rnaseq/batch_gene_reports.R",
    genes=ancient("reports/rnaseq/genes_to_report.txt")
  params:
    outdir="reports/rnaseq/gene_reports"
  #output:
  #  directory("reports/rnaseq/gene_reports")
  shell:
    "Rscript {input.script}"
    " --genes {input.genes}"
    " --outdir {params.outdir}"

#### BCR and antibody isotype analyses ####

# Download TRUST to projdir/bin
# Can be accessed via git clone https://github.com/liulab-dfci/TRUST4.git
# TRUST sometimes fails when invoked directly with Snakemake
# In this case, try src/rnaseq/trust4.sh

rule trust:
  input:
    fq="data/rnaseq/fastq/{sample}.fastq.gz",
    bcrtcrfa=ancient("bin/TRUST4/hg38_bcrtcr.fa"),
    imgt=ancient("bin/TRUST4/human_IMGT+C.fa"),
    trust=ancient("bin/TRUST4/run-trust4")
  output:
    report="data/rnaseq/TRUST/{sample}.fastq.gz_report.tsv"
  params:
    threads=16,
    outdir="data/rnaseq/TRUST/"
  shell:
    """
    mkdir -p {params.outdir} &&
    {input.trust} -u {input.fq} -t {params.threads} -f {input.bcrtcrfa} --ref {input.imgt} -o {params.outdir}{input.fq}
    """

rule trust_report:
  input:
    fqdata="data/rnaseq/metadata/01_rnaMeta.Rds",
    sampledata="data/rnaseq/metadata/05_sample_annot_filtered.csv",
    survdata="data/rnaseq/interim/04_survdata.Rds",
    preexcluded_samples="data/external/pre_excluded_samples.csv",
    discarded="data/rnaseq/metadata/02_discarded_samples.csv",
    ihc_outliers="data/rnaseq/metadata/03_removed_pam50_outliers.csv",
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    cp="data/rnaseq/interim/color_palettes.Rds",
    rmd="reports/rnaseq/15_BCR_clonality.Rmd",
    script="src/utils/rmarkdown.R"
  params:
    trustdir="data/rnaseq/TRUST"
  output:
    alltrust="data/rnaseq/processed/15_alltrust.Rds",
    trustdata="data/rnaseq/processed/15_trustdata.Rds",
    trustexcel="results/rnaseq/TRUST/TRUST_results.xlsx",
    km_cpk="figures/supfigs/Supfig21_clonality_km_OS_DRS.pdf",
    html="reports/rnaseq/15_BCR_clonality.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --trustdir {params.trustdir}"
    " --fqdata {input.fqdata}"
    " --sampledata {input.sampledata}"
    " --survdata {input.survdata}"
    " --preexcluded_samples {input.preexcluded_samples}"
    " --discarded {input.discarded}"
    " --ihc_outliers {input.ihc_outliers}"
    " --dds {input.dds}"
    " --cp {input.cp}"
    " --km_cpk {output.km_cpk}"
    
rule antibody_isotypes:
  input:
    script="src/utils/rmarkdown.R",
    rmd="reports/rnaseq/15b_antibody_isotypes.Rmd",
    bx_annot=ancient("data/rnaseq/processed/bx_annot.Rds"),
    dds=ancient("data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds")
  output:
    html="reports/rnaseq/15b_antibody_isotypes.html",
    igNorm="data/processed/15b_igNorm.Rds",
    igSurv="data/processed/15b_igSurv.Rds"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"  
    " --dds {input.dds}"
    " --bx_annot {input.bx_annot}"

#### Shiny App for PPBC gene visualization ####

# Set up the data to be loaded within the shiny app
rule app_setup:
  input:
    diffex=expand("results/rnaseq/diffex/{res}.xlsx", 
      res=["08_one_vs_rest_allgenes", "07_pairwise_comparisons_allgenes", "06_LRT_allgenes"]),
    bx_annot="data/rnaseq/processed/bx_annot.Rds",
    coxres="results/rnaseq/survival/12_cox_allgenes.xlsx",
    duplicates="src/rnaseq/summarize_duplicate_ids.R",
    coxdata="data/rnaseq/processed/12_coxdata.Rds"
  output:
    "shinyApp/VisualizePPBCgene/data/app_gx_annot.Rds",
    "shinyApp/VisualizePPBCgene/data/12_cox_allgenes.xlsx",
    "shinyApp/VisualizePPBCgene/data/app_diffex_res_list.Rds",
    "shinyApp/VisualizePPBCgene/data/app_survival_sample_data.Rds",
    "shinyApp/VisualizePPBCgene/data/app_ensembl_tmmnorm_genesxsample.Rds",
    "shinyApp/VisualizePPBCgene/data/app_symbol_tmmnorm_genesxsample.Rds"
  shell:
    """
    Rscript shinyApp/VisualizePPBCgene/data_setup.R
    """ 

