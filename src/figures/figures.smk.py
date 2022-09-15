rule aggregate_results:
  input:
    patientdata="data/external/patient_data.tsv",
    dds="data/rnaseq/interim/05_dds_PAM50_batch.Rds",
    ciber_input="results/rnaseq/cibersortX/CIBERSORTx_Job6_Results.csv",
    igSurv="data/processed/15b_igSurv.Rds",
    cluster_results="results/rnaseq/clustering/11_inv_clusters.xlsx", #for column clusters
    read_patient_data="src/utils/read_patient_data.R",
    quantile_categories="src/figures/quantile_categories.R",
    rmd="reports/figures/00_aggregate_results.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    html="reports/figures/00_aggregate_results.html",
    figuredata="data/figures/00_figuredata.Rds"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.html}"
    " --patientdata {input.patientdata}"
    " --dds {input.dds}"
    " --ciber_input {input.ciber_input}"
    " --cluster_results {input.cluster_results}"
    " --igSurv {input.igSurv}"
    " --read_patient_data {input.read_patient_data}"
    " --quantile_categories {input.quantile_categories}"
    
rule ig_clusters:
  input:
    cluster_results="results/rnaseq/clustering/11_inv_clusters.xlsx", #for row clusters
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    bx_annot="data/rnaseq/processed/bx_annot.Rds",
    figuredata="data/figures/00_figuredata.Rds",
    rmd="reports/figures/01_IG_clusters.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    report="reports/figures/01_IG_clusters.pdf",
    fig1b="figures/Fig1b_ig_boxplot.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --cluster_results {input.cluster_results}"
    " --dds {input.dds}"
    " --bx_annot {input.bx_annot}"
    " --figuredata {input.figuredata}"
    " --fig1b {output.fig1b}"

rule spatial_figures:
  input:
    figuredata="data/figures/00_figuredata.Rds",
    sampledata = "data/external/sample_data.tsv",
    density="data/vectra/processed/density_ppbc.Rds",
    spatstat_l="data/vectra/processed/11_spatstat_l.Rds",
    spatstat_lcross_immune="data/vectra/processed/11_spatstat_lcross_immune.Rds",
    spatstat_lcross_panck="data/vectra/processed/11_spatstat_lcross_panck.Rds",
    rmd="reports/figures/02_spatial_figures.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    report="reports/figures/02_spatial_figures.pdf",
    fig2b="figures/Fig2b_CD20_density_boxplot.pdf",
    fig2c="figures/Fig2c_CD20_density_km.pdf",
    rt_fig2c="figures/Fig2c_risktable.csv",
    fig2e="figures/Fig2e_CD20_l_boxplot.pdf",
    fig2f="figures/Fig2f_Cd20_l_km.pdf",
    rt_fig2f="figures/Fig2f_risktable.csv"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --figuredata {input.figuredata}"
    " --sampledata {input.sampledata}"
    " --density {input.density}"
    " --spatstat_l {input.spatstat_l}"
    " --spatstat_lcross_immune {input.spatstat_lcross_immune}"
    " --spatstat_lcross_panck {input.spatstat_lcross_panck}"
    " --fig2b {output.fig2b}"
    " --fig2c {output.fig2c}"
    " --rt_fig2c {output.rt_fig2c}"
    " --fig2e {output.fig2e}"
    " --fig2f {output.fig2f}"
    " --rt_fig2f {output.rt_fig2f}"
    
#vsd="data/rnaseq/interim/08_vsd_ovr.Rds",
#cluster_survival="data/rnaseq/interim/11_ig_survdata.Rds",
