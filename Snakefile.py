#### General use ####

# Type "snakemake" with no arguments into console to generate the files listed under rule all
# Alternatively, specify a rule or output file

# Example, performs a dry run of the gene_reports rule with 8 cores
# snakemake -np gene_reports -s Snakefile.py -j 8

# Useful arguments
# -j designates thread number
# -s specifies a Snakefile path other than the default ("Snakefile" in working dir)
# -n for dry run, no rule execution
# -p print shell commands of rules that will be executed, useful together with -n

# The Snakefile is given a .py extension for automatic syntax highlighting in IDEs
# Either use snakemake -s Snakefile.py, or create a symlink (Snakefile -> Snakefile.py)

# To get a graphical representation of the workflow:
# "snakemake --dag | dot -Tsvg > dag.svg"

# Subworkflows can be included as additional Snakefiles that are topic specific
# These are designated with a .smk.py suffix

#### Conda and Snakemake ####

# Create a new conda environment for a specific software version with: 
#  conda create -n myenv python=3.9

# Save environment to text with: 
#  conda activate myenv; conda env export > envs/myenv.yml

# For a condensed version:
#  conda env export --from-history > envs/myenv_brief.yml

# To update an existing environment file: 
#  conda env update --prefix ./envs --file myenv.yml  --prune

# You can optionally remove the environment from conda once the yml file has been created, 
# Snakemake will create it anew when rules are run

# To activate conda environment from a config file, use:
#  conda env create -f envs/myenv.yml

# To use a conda envinroment with snakemake: 
#  snakemake --use-conda

# Note: the environment file path is relative to the (sub-)Snakefile

# To create the environments without running any rules: 
#  snakemake --use-conda --create-envs-only

from pathlib import Path

configfile: "config.yaml"

include: "src/rnaseq/rnaseq.smk.py"
include: "src/spatial/spatial.smk.py"
include: "src/figures/figures.smk.py"

rule all:
  input:
    # Figures
    "data/figures/00_figuredata.Rds",
    "figures/Fig2b_CD20_density_boxplot.pdf",
    "figures/Fig2c_CD20_density_km_OS.pdf",
    "figures/Fig2e_CD20_l_boxplot.pdf",
    "figures/Fig2f_Cd20_l_km_OS.pdf",
    "figures/Fig2h_TAPC_PPBC_cor.pdf",
    "figures/Fig2i_TAPC_KM_OS.pdf",
    "figures/Fig2k_CD38_boxplot.pdf",
    "figures/Fig2l_CD38_KM_OS.pdf",
    "figures/Fig3a_boxplot_isotypes.pdf",
    "figures/Fig3bc_km_isotypes.pdf",
    "figures/Fig3d_cor_Ig_TAPC.pdf",
    "figures/Fig4a_TIL_boxplot.pdf",
    "figures/Fig4b_cibersort_CD8_boxplot.pdf",
    "figures/Fig4c_cibersort_CD8_km.pdf",
    "figures/Fig4d_CD8_density_boxplot_km_OS.pdf",
    "figures/Fig4e_CD8_l_boxplot_km_OS.pdf",
    "figures/Fig4f_CD8_PanCk_lcross_boxplot_km_OS.pdf",
    "figures/Fig4g_CD8_CD4_lcross_boxplot_km_OS.pdf",
    "figures/Fig4h_CD20_CD8_lcross_boxplot_km_OS.pdf",
    "figures/Fig4i_CD20_CD4_lcross_boxplot_km_OS.pdf",
    # Sup Figs
    "figures/supfigs/Supfig2_studygroup_km.pdf",
    "figures/supfigs/Supfig3_pca_rnaseq.pdf",
    "figures/supfigs/Supfig4a_volc_inv_rest.pdf",
    "figures/supfigs/Supfig4b_volc_inv_nonprbc.pdf",
    "figures/supfigs/Supfig4c_volc_inv_prbc.pdf",
    "figures/supfigs/Supfig4d_volc_inv_lac.pdf",
    "figures/supfigs/Supfig5a_volc_involution_duration.pdf",
    "figures/supfigs/Supfig5b_volc_breastfeeding_duration.pdf",
    "figures/supfigs/Supfig6d_volc_ppbcpw_vs_rest_subgroup_Basal.pdf",
    "figures/supfigs/Supfig6c_volc_ppbcpw_vs_rest_subgroup_Her2.pdf",
    "figures/supfigs/Supfig6a_volc_ppbcpw_vs_rest_subgroup_LumA.pdf",
    "figures/supfigs/Supfig6b_volc_ppbcpw_vs_rest_subgroup_LumB.pdf",
    "figures/supfigs/Supfig7a_milk_vs_IG_genes.pdf",
    "figures/supfigs/Supfig7b_cor_ig_milk_all.pdf",
    "figures/supfigs/Supfig7c_cor_ig_milk_inv.pdf",
    "figures/supfigs/Supfig8_igSig_km_forest.pdf",
    "figures/supfigs/Supfig9_IG_cluster_PPBC_KM_DRS.pdf",
    "figures/supfigs/Supfig10a_IG_cluster_PAM50_KM_OS.pdf",
    "figures/supfigs/Supfig10b_IG_cluster_PAM50_KM_DRS.pdf",
    "figures/supfigs/Supfig11a_uni_genewise_drs_heatmap.pdf",
    #The panels in Supfig11 are all derived from the reports generated by rule figure_genes
    "figures/supfigs/SupFig12a_cibersort_relative_heatmap.pdf",
    "figures/supfigs/Supfig12b_boxplot_cibersort_B_subtypes.pdf",
    "figures/supfigs/Supfig12c_kaplan_cibersort_plasmaB_OS_DRS.pdf",
    "figures/supfigs/Supfig12d_kaplan_cibersort_memoryB_OS_DRS.pdf",
    "figures/supfigs/Supfig13a_CD20_density_km_DRS.pdf",
    "figures/supfigs/Supfig13b_CD20_l_km_DRS.pdf",
    "figures/supfigs/Supfig13c_TAPC_KM_DRS.pdf",
    "figures/supfigs/Supfig13d_CD38_KM_DRS.pdf",
    "figures/supfigs/Supfig14a_CD20_Panck_lcross_boxplot.pdf",
    "figures/supfigs/Supfig14b_CD20_Panck_lcross_km_OS.pdf",
    "figures/supfigs/Supfig14c_CD20_Panck_lcross_km_DRS.pdf",
    "figures/supfigs/Supfig16_TAPC_CD38_cor.pdf",
    "figures/supfigs/Supfig17_cor_IgA_milk.pdf",
    "figures/supfigs/Supfig18_km_isotypes_drs.pdf",
    "figures/supfigs/Supfig19_TAPC_isotype_cor.pdf",
    "figures/supfigs/Supfig20_cor_IgA_IgG.pdf",
    "figures/supfigs/Supfig21_clonality_km_OS_DRS.pdf",
    "figures/supfigs/Supfig23A_cibersort_CD8_km_DRS.pdf",
    "figures/supfigs/Supfig23b-g_misc_spatial_km_drs.pdf",
    "figures/supfigs/Supfig27_boxplot_cibersort_other_cells.pdf",
    "figures/supfigs/Supfig28_spatial_boxplot_other_celltypes.pdf"

# Utility for converting Excel metadata to text
# Text metadata can be tracked via git (if it's not too large)
rule excel_to_tsv:
  input:
    script="src/utils/excel_to_tsv.R",
    lib="src/utils/parse_args.R",
    metadata="data/external/PPBC_metadata_20220811.xlsx"
  params:
    outDir="data/external"
  output:
    sampledata = "data/external/sample_data.tsv",
    patientdata = "data/external/patient_data.tsv",
    codebook = "data/external/codebook.tsv"
  shell:
    "Rscript {input.script} "
    "--metadata {input.metadata} "
    "--outDir {params.outDir}"
