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
    report="reports/figures/IG_clusters_km.pdf",
    fig1b="figures/Fig1b_ig_boxplot.pdf",
    supfig2="figures/supfigs/Supfig2_studygroup_km.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --cluster_results {input.cluster_results}"
    " --dds {input.dds}"
    " --bx_annot {input.bx_annot}"
    " --figuredata {input.figuredata}"
    " --fig1b {output.fig1b}"
    " --supfig2 {output.supfig2}"

rule spatial_figures:
  input:
    figuredata="data/figures/00_figuredata.Rds",
    sampledata = "data/external/sample_data.tsv",
    density="data/vectra/processed/density_ppbc.Rds",
    spatstat_l="data/vectra/processed/11_spatstat_l.Rds",
    spatstat_lcross_immune="data/vectra/processed/11_spatstat_lcross_immune.Rds",
    spatstat_lcross_panck="data/vectra/processed/11_spatstat_lcross_panck.Rds",
    rmd="reports/figures/spatial_figures.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    report="reports/figures/spatial_figures.pdf",
    fig2b="figures/Fig2b_CD20_density_boxplot.pdf",
    fig2c="figures/Fig2c_CD20_density_km.pdf",
    rt_fig2c="figures/Fig2c_risktable.csv",
    fig2e="figures/Fig2e_CD20_l_boxplot.pdf",
    fig2f="figures/Fig2f_Cd20_l_km.pdf",
    rt_fig2f="figures/Fig2f_risktable.csv",
    fig4d="figures/Fig4d_CD8_density_boxplot_km.pdf",
    rt_fig4d="figures/Fig4d_risktable.csv",
    fig4e="figures/Fig4e_CD8_l_boxplot_km.pdf",
    rt_fig4e="figures/Fig4e_risktable.csv",
    fig4f="figures/Fig4f_CD8_PanCk_lcross_boxplot_km.pdf",
    rt_fig4f="figures/Fig4f_risktable.csv",
    fig4g="figures/Fig4g_CD8_CD4_lcross_boxplot_km.pdf",
    rt_fig4g="figures/Fig4g_risktable.csv",
    fig4h="figures/Fig4h_CD20_CD8_lcross_boxplot_km.pdf",
    rt_fig4h="figures/Fig4h_risktable.csv",
    fig4i="figures/Fig4i_CD20_CD4_lcross_boxplot_km.pdf",
    rt_fig4i="figures/Fig4i_risktable.csv"
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
    " --fig4d {output.fig4d}"
    " --rt_fig4d {output.rt_fig4d}"
    " --fig4e {output.fig4e}"
    " --rt_fig4e {output.rt_fig4e}"
    " --fig4f {output.fig4f}"
    " --rt_fig4f {output.rt_fig4f}"
    " --fig4g {output.fig4g}"
    " --rt_fig4g {output.rt_fig4g}"
    " --fig4h {output.fig4h}"
    " --rt_fig4h {output.rt_fig4h}"
    " --fig4i {output.fig4i}"
    " --rt_fig4i {output.rt_fig4i}"

rule stainings:
  input:
    figuredata="data/figures/00_figuredata.Rds",
    rmd="reports/figures/CD38_TAPC_TILs.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    report="reports/figures/CD38_TAPC_TILs.pdf",
    fig2h="figures/Fig2h_TAPC_PPBC_cor.pdf",
    fig2i="figures/Fig2i_TAPC_KM.pdf",
    rt_fig2i="figures/Fig2i_risktable.csv",
    fig2k="figures/Fig2k_CD38_boxplot.pdf",
    fig2l="figures/Fig2l_CD38_KM.pdf",
    rt_fig2l="figures/Fig2l_risktable.csv",
    fig4a="figures/Fig4a_TIL_boxplot.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --figuredata {input.figuredata}"
    " --fig2h {output.fig2h}"
    " --fig2i {output.fig2i}"
    " --rt_fig2i {output.rt_fig2i}"
    " --fig4a {output.fig4a}"

rule isotypes:
  input:
    figuredata="data/figures/00_figuredata.Rds",
    rmd="reports/figures/isotypes.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    fig3a="figures/Fig3a_boxplot_isotypes.pdf",
    fig3bc="figures/Fig3bc_km_isotypes.pdf",
    rt_fig3bc="figures/Fig3bc_risktable.csv",
    fig3d="figures/Fig3d_cor_Ig_TAPC.pdf",
    report="reports/figures/isotypes.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --figuredata {input.figuredata}"
    " --fig3a {output.fig3a}"
    " --fig3bc {output.fig3bc}"
    " --rt_fig3bc {output.rt_fig3bc}"
    " --fig3d {output.fig3d}"
    
rule ciberfigs:
  input:
    figuredata="data/figures/00_figuredata.Rds",
    rmd="reports/figures/ciberfigs.Rmd",
    script="src/utils/rmarkdown.R"
  output:
    fig4b="figures/Fig4b_cibersort_CD8_boxplot.pdf",
    fig4c="figures/Fig4c_cibersort_CD8_km.pdf",
    rt_fig4c="figures/Fig4c_risktable.csv",
    report="reports/figures/ciberfigs.pdf"
  shell:
    "Rscript {input.script} {input.rmd} $PWD/{output.report}"
    " --figuredata {input.figuredata}"
    " --fig4b {output.fig4b}"
    " --fig4c {output.fig4c}"
    " --rt_fig4c {output.rt_fig4c}"
    
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
    cor_ig_milk_inv="results/rnaseq/diffex/cor_ig_milk_inv.pdf"
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
    cor_ig_milk_inv="figures/supfigs/Supfig7c_cor_ig_milk_inv.pdf"
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
    """
    
    
#vsd="data/rnaseq/interim/08_vsd_ovr.Rds",
#cluster_survival="data/rnaseq/interim/11_ig_survdata.Rds",
