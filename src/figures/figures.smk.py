rule aggregate_results:
  input:
    patientdata="data/external/patient_data.tsv",
    dds="data/rnaseq/interim/05_dds_PAM50_batch.Rds",
    ciber_input="results/rnaseq/cibersortX/CIBERSORTx_Job6_Results.csv",
    igSurv="data/processed/15b_igSurv.Rds",
    cluster_results="results/rnaseq/clustering/11_inv_clusters.xlsx", #for column clusters
    read_patient_data="src/utils/read_patient_data.R",
    quantile_categories="src/figures/quantile_categories.R",
    rmd="reports/figures/aggregate_results.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    html="reports/figures/aggregate_results.html",
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
    
rule ig_clusters_km:
  input:
    cluster_results="results/rnaseq/clustering/11_inv_clusters.xlsx", #for row clusters
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    bx_annot="data/rnaseq/processed/bx_annot.Rds",
    figuredata="data/figures/00_figuredata.Rds",
    rmd="reports/figures/IG_clusters_km.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    report="reports/figures/IG_clusters_km.html",
    fig1b="figures/Fig1b_ig_boxplot.pdf",
    supfig2="figures/supfigs/Supfig2_studygroup_km.pdf",
    supfig8="figures/supfigs/Supfig8_igSig_km_forest.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --cluster_results {input.cluster_results}"
    " --dds {input.dds}"
    " --bx_annot {input.bx_annot}"
    " --figuredata {input.figuredata}"
    " --fig1b {output.fig1b}"
    " --supfig2 {output.supfig2}"
    " --supfig8 {output.supfig8}"

rule spatial_figures:
  input:
    figuredata="data/figures/00_figuredata.Rds",
    sampledata = "data/external/sample_data.tsv",
    density="data/vectra/processed/density_ppbc.Rds",
    spatstat_l="data/vectra/processed/11_spatstat_l.Rds",
    spatstat_lcross_immune="data/vectra/processed/11_spatstat_lcross_immune.Rds",
    spatstat_lcross_panck="data/vectra/processed/11_spatstat_lcross_panck.Rds",
    rt="src/figures/faceted_risktable.R",
    rmd="reports/figures/spatial_figures.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    report="reports/figures/spatial_figures.html",
    boxplot_CD20_density="figures/Fig2b_CD20_density_boxplot.pdf",
    km_CD20_density_OS="figures/Fig2c_CD20_density_km_OS.pdf",
    boxplot_CD20_lstat="figures/Fig2e_CD20_l_boxplot.pdf",
    km_CD20_lstat_OS="figures/Fig2f_Cd20_l_km_OS.pdf",
    boxplot_CD8_density_km_OS="figures/Fig4d_CD8_density_boxplot_km_OS.pdf",
    boxplot_CD8_l_km_OS="figures/Fig4e_CD8_l_boxplot_km_OS.pdf",
    boxplot_CD8_PancK_lcross_km_OS="figures/Fig4f_CD8_PanCk_lcross_boxplot_km_OS.pdf",
    boxplot_CD8_CD4_lcross_km_OS="figures/Fig4g_CD8_CD4_lcross_boxplot_km_OS.pdf",
    boxplot_CD20_CD8_lcross_km_OS="figures/Fig4h_CD20_CD8_lcross_boxplot_km_OS.pdf",
    boxplot_CD20_CD4_lcross_km_OS="figures/Fig4i_CD20_CD4_lcross_boxplot_km_OS.pdf",
    km_CD20_density_DRS="figures/supfigs/Supfig13a_CD20_density_km_DRS.pdf",
    km_CD20_l_DRS="figures/supfigs/Supfig13b_CD20_l_km_DRS.pdf",
    boxplot_CD20_Panck_lcross="figures/supfigs/Supfig14a_CD20_Panck_lcross_boxplot.pdf",
    km_CD20_Panck_lcross_km_OS="figures/supfigs/Supfig14b_CD20_Panck_lcross_km_OS.pdf",
    km_CD20_Panck_lcross_km_DRS="figures/supfigs/Supfig14c_CD20_Panck_lcross_km_DRS.pdf",
    misc_spatial_drs="figures/supfigs/Supfig23b-g_misc_spatial_km_drs.pdf",
    boxplot_others="figures/supfigs/Supfig28_spatial_boxplot_other_celltypes.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --figuredata {input.figuredata}"
    " --sampledata {input.sampledata}"
    " --density {input.density}"
    " --spatstat_l {input.spatstat_l}"
    " --spatstat_lcross_immune {input.spatstat_lcross_immune}"
    " --spatstat_lcross_panck {input.spatstat_lcross_panck}"
    " --boxplot_CD20_density {output.boxplot_CD20_density}"
    " --km_CD20_density_OS {output.km_CD20_density_OS}"
    " --boxplot_CD20_lstat {output.boxplot_CD20_lstat}"
    " --km_CD20_lstat_OS {output.km_CD20_lstat_OS}"
    " --boxplot_CD8_density_km_OS {output.boxplot_CD8_density_km_OS}"
    " --boxplot_CD8_l_km_OS {output.boxplot_CD8_l_km_OS}"
    " --boxplot_CD8_PancK_lcross_km_OS {output.boxplot_CD8_PancK_lcross_km_OS}"
    " --boxplot_CD8_CD4_lcross_km_OS {output.boxplot_CD8_CD4_lcross_km_OS}"
    " --boxplot_CD20_CD8_lcross_km_OS {output.boxplot_CD20_CD8_lcross_km_OS}"
    " --boxplot_CD20_CD4_lcross_km_OS {output.boxplot_CD20_CD4_lcross_km_OS}"
    " --km_CD20_density_DRS {output.km_CD20_density_DRS}"
    " --km_CD20_l_DRS {output.km_CD20_l_DRS}"
    " --boxplot_CD20_Panck_lcross {output.boxplot_CD20_Panck_lcross}"
    " --km_CD20_Panck_lcross_km_OS {output.km_CD20_Panck_lcross_km_OS}"
    " --km_CD20_Panck_lcross_km_DRS {output.km_CD20_Panck_lcross_km_DRS}"
    " --misc_spatial_drs {output.misc_spatial_drs}"
    " --boxplot_others {output.boxplot_others}"

rule stainings:
  input:
    figuredata="data/figures/00_figuredata.Rds",
    rmd="reports/figures/CD38_TAPC_TILs.Rmd",
    rt="src/figures/faceted_risktable.R",
    study_km="src/figures/study_km.R",
    script="src/utils/rmarkdown.R"
  output:
    report="reports/figures/CD38_TAPC_TILs.html",
    TAPC_PPBC_cor="figures/Fig2h_TAPC_PPBC_cor.pdf",
    km_TAPC_OS="figures/Fig2i_TAPC_KM_OS.pdf",
    boxplot_CD38="figures/Fig2k_CD38_boxplot.pdf",
    km_CD38_OS="figures/Fig2l_CD38_KM_OS.pdf",
    boxplot_TIL="figures/Fig4a_TIL_boxplot.pdf",
    km_TAPC_DRS="figures/supfigs/Supfig13c_TAPC_KM_DRS.pdf",
    km_CD38_DRS="figures/supfigs/Supfig13d_CD38_KM_DRS.pdf",
    TAPC_CD38_cor="figures/supfigs/Supfig16_TAPC_CD38_cor.pdf",
    TAPC_iso_cor="figures/supfigs/Supfig19_TAPC_isotype_cor.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --figuredata {input.figuredata}"
    " --TAPC_PPBC_cor {output.TAPC_PPBC_cor}"
    " --km_TAPC_OS {output.km_TAPC_OS}"
    " --boxplot_CD38 {output.boxplot_CD38}"
    " --km_CD38_OS {output.km_CD38_OS}"
    " --boxplot_TIL {output.boxplot_TIL}"
    " --km_TAPC_DRS {output.km_TAPC_DRS}"
    " --km_CD38_DRS {output.km_CD38_DRS}"
    " --TAPC_CD38_cor {output.TAPC_CD38_cor}"
    " --TAPC_iso_cor {output.TAPC_iso_cor}"

rule isotypes:
  input:
    figuredata="data/figures/00_figuredata.Rds",
    normat="data/rnaseq/processed/08_sfnorm_genesymbol_countmatrix.Rds",
    ig_milk_genes="data/rnaseq/processed/08_ig_and_milk_genes.Rds",
    rmd="reports/figures/isotypes.Rmd",
    rt="src/figures/faceted_risktable.R",
    script="src/utils/rmarkdown.R"
  output:
    boxplot_abs="figures/Fig3a_boxplot_isotypes.pdf",
    km_ig_os="figures/Fig3bc_km_isotypes.pdf",
    km_ig_drs="figures/supfigs/Supfig18_km_isotypes_drs.pdf",
    cor_IG_TAPC="figures/Fig3d_cor_Ig_TAPC.pdf",
    cor_IgA_milk="figures/supfigs/Supfig17_cor_IgA_milk.pdf",
    cor_IgA_IgG="figures/supfigs/Supfig20_cor_IgA_IgG.pdf",
    report="reports/figures/isotypes.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --figuredata {input.figuredata}"
    " --normat {input.normat}"
    " --ig_milk_genes {input.ig_milk_genes}"
    " --rt {input.rt}"
    " --boxplot_abs {output.boxplot_abs}"
    " --km_ig_os {output.km_ig_os}"
    " --km_ig_drs {output.km_ig_drs}"
    " --cor_IG_TAPC {output.cor_IG_TAPC}"
    " --cor_IgA_milk {output.cor_IgA_milk}"
    " --cor_IgA_IgG {output.cor_IgA_IgG}"
    
rule ciberfigs:
  input:
    figuredata="data/figures/00_figuredata.Rds",
    rmd="reports/figures/ciberfigs.Rmd",
    rt="src/figures/faceted_risktable.R",
    script="src/utils/rmarkdown.R"
  output:
    boxplot_cd8="figures/Fig4b_cibersort_CD8_boxplot.pdf",
    km_cd8="figures/Fig4c_cibersort_CD8_km.pdf",
    boxplot_b_subtypes="figures/supfigs/Supfig12b_boxplot_cibersort_B_subtypes.pdf",
    km_plasmaB="figures/supfigs/Supfig12c_kaplan_cibersort_plasmaB_OS_DRS.pdf",
    km_memoryB="figures/supfigs/Supfig12d_kaplan_cibersort_memoryB_OS_DRS.pdf",
    km_cd8_DRS="figures/supfigs/Supfig23A_cibersort_CD8_km_DRS.pdf",
    bp_other="figures/supfigs/Supfig27_boxplot_cibersort_other_cells.pdf",
    report="reports/figures/ciberfigs.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --figuredata {input.figuredata}"
    " --boxplot_cd8 {output.boxplot_cd8}"
    " --km_cd8 {output.km_cd8}"
    " --boxplot_b_subtypes {output.boxplot_b_subtypes}"
    " --km_plasmaB {output.km_plasmaB}"
    " --km_memoryB {output.km_memoryB}"
    " --km_cd8_DRS {output.km_cd8_DRS}"
    " --bp_other {output.bp_other}"
    
rule figure_genes:
  input:
    aggdata="data/rnaseq/processed/16_gene_report_environment.RData",
    bx_annot="data/rnaseq/processed/bx_annot.Rds",
    dds="data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds",
    rmd="reports/rnaseq/17_gene_report_template.Rmd",
    script="src/rnaseq/batch_gene_reports.R",
    genes="reports/figures/figure_genes.txt",
  params:
    outdir="figures/figure_genes"
  # output:
  #   directory("reports/rnaseq/figure_genes")
  shell:
    "Rscript {input.script}"
    " --genes {input.genes}"
    " --outdir {params.outdir}"
    
rule reviewer_requests:
  input:
    figuredata="data/figures/00_figuredata.Rds",
    rmd="reports/figures/reviewer_requests.Rmd",
    rt="src/figures/faceted_risktable.R",
    script="src/utils/rmarkdown.R"
  output:
    report="reports/figures/reviewer_requests.html",
    reviewer_requests="figures/reviewer_requests.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --figuredata {input.figuredata}"
    " --reviewer_requests {output.reviewer_requests}"

rule copy_figs:
  input:
    rnaseq_pcas="results/rnaseq/dimensionality/pca_rnaseq.pdf",
    vp_inv_rest="results/rnaseq/diffex/volc_inv_rest.pdf",
    vp_inv_nonprbc="results/rnaseq/diffex/volc_inv_nonprbc.pdf",
    vp_inv_prbc="results/rnaseq/diffex/volc_inv_prbc.pdf",
    vp_inv_lac="results/rnaseq/diffex/volc_inv_nonprbc.pdf",
    vp_duration_inv="results/rnaseq/diffex/09_volcano_involution_duration.pdf",
    vp_duration_bf="results/rnaseq/diffex/09b_volcano_breastfeeding_duration.pdf",
    vp_sub_invrest_basal="results/rnaseq/diffex/14_volc_ppbcpw_vs_rest_subgroup_Basal.pdf",
    vp_sub_invrest_her2="results/rnaseq/diffex/14_volc_ppbcpw_vs_rest_subgroup_Her2.pdf",
    vp_sub_invrest_lumA="results/rnaseq/diffex/14_volc_ppbcpw_vs_rest_subgroup_LumA.pdf",
    vp_sub_invrest_lumB="results/rnaseq/diffex/14_volc_ppbcpw_vs_rest_subgroup_LumB.pdf",
    milk_heatmap="results/rnaseq/diffex/milk_vs_IG_genes.pdf",
    cor_ig_milk_all="results/rnaseq/diffex/cor_ig_milk_all.pdf",
    cor_ig_milk_inv="results/rnaseq/diffex/cor_ig_milk_inv.pdf",
    uni_genewise_drs_heatmap="results/rnaseq/survival/12_uni_genewise_drs_heatmap.pdf"
  output:
    rnaseq_pcas="figures/supfigs/Supfig3_pca_rnaseq.pdf",
    vp_inv_rest="figures/supfigs/Supfig4a_volc_inv_rest.pdf",
    vp_inv_nonprbc="figures/supfigs/Supfig4b_volc_inv_nonprbc.pdf",
    vp_inv_prbc="figures/supfigs/Supfig4c_volc_inv_prbc.pdf",
    vp_inv_lac="figures/supfigs/Supfig4d_volc_inv_lac.pdf",
    vp_duration_inv="figures/supfigs/Supfig5a_volc_involution_duration.pdf",
    vp_duration_bf="figures/supfigs/Supfig5b_volc_breastfeeding_duration.pdf",
    vp_sub_invrest_basal="figures/supfigs/Supfig6d_volc_ppbcpw_vs_rest_subgroup_Basal.pdf",
    vp_sub_invrest_her2="figures/supfigs/Supfig6c_volc_ppbcpw_vs_rest_subgroup_Her2.pdf",
    vp_sub_invrest_lumA="figures/supfigs/Supfig6a_volc_ppbcpw_vs_rest_subgroup_LumA.pdf",
    vp_sub_invrest_lumB="figures/supfigs/Supfig6b_volc_ppbcpw_vs_rest_subgroup_LumB.pdf",
    milk_heatmap="figures/supfigs/Supfig7a_milk_vs_IG_genes.pdf",
    cor_ig_milk_all="figures/supfigs/Supfig7b_cor_ig_milk_all.pdf",
    cor_ig_milk_inv="figures/supfigs/Supfig7c_cor_ig_milk_inv.pdf",
    uni_genewise_drs_heatmap="figures/supfigs/Supfig11a_uni_genewise_drs_heatmap.pdf"
  shell:
    """
    cp -v {input.rnaseq_pcas} {output.rnaseq_pcas}
    cp -v {input.vp_inv_rest} {output.vp_inv_rest}
    cp -v {input.vp_inv_nonprbc} {output.vp_inv_nonprbc}
    cp -v {input.vp_inv_prbc} {output.vp_inv_prbc}
    cp -v {input.vp_inv_lac} {output.vp_inv_lac}
    cp -v {input.vp_duration_inv} {output.vp_duration_inv}
    cp -v {input.vp_duration_bf} {output.vp_duration_bf}
    cp -v {input.vp_sub_invrest_basal} {output.vp_sub_invrest_basal}
    cp -v {input.vp_sub_invrest_her2} {output.vp_sub_invrest_her2}
    cp -v {input.vp_sub_invrest_lumA} {output.vp_sub_invrest_lumA}
    cp -v {input.vp_sub_invrest_lumB} {output.vp_sub_invrest_lumB}
    cp -v {input.milk_heatmap} {output.milk_heatmap}
    cp -v {input.cor_ig_milk_all} {output.cor_ig_milk_all}
    cp -v {input.cor_ig_milk_inv} {output.cor_ig_milk_inv}
    cp -v {input.uni_genewise_drs_heatmap} {output.uni_genewise_drs_heatmap}
    """
    
rule renumber_figs:
  input:
    deg_hm="figures/Fig1a_DEG_heatmap.pdf",
    deg_hm2="figures/Fig1b_DEG_heatmap_kmeans.pdf",
    igboxplot="figures/Fig1b_ig_boxplot.pdf",
    igbarplot="figures/Fig1e_IG_cluster_barplot.pdf",
    igclustkm="figures/Fig1f_IG_cluster_PPBC_KM_OS.pdf",
    rmd="reports/figures/figure_numbering.Rmd",
    fileLoc="reports/figures/figure_numbering.csv",
    script="src/utils/rmarkdown.R"
  output:
    deg_hm="figures/Fig2a_DEG_heatmap.pdf",
    deg_hm2="figures/Fig2b_DEG_heatmap_kmeans.pdf",
    igboxplot="figures/Fig2c_ig_boxplot.pdf",
    igbarplot="figures/Fig2e_IG_cluster_barplot.pdf",
    igclustkm="figures/Fig2f_IG_cluster_PPBC_KM_OS.pdf",
    report="reports/figures/figure_numbering.html"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --fileLoc {input.fileLoc}"
