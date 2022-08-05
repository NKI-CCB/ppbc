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
#To activate conda environment, use: conda env create -f envs/environment.yml
#If environment.yml has been updated, update the environment with: 
#conda env update --prefix ./envs --file environment.yml  --prune
#Save environment to text with: conda env export > envs/environment.yml
#For condensed version: conda env export --from-history > envs/env_brief.yml
#To use a conda envinroment with snakemake: snakemake -n --use-conda
#To create the environments without running any rules: 
#snakemake -n --use-conda --create-envs-only

from pathlib import Path

configfile: "config.yaml"

include: "src/rnaseq/rnaseq.smk.py"
include: "src/spatial/spatial.smk.py"

rule all:
  input:
    # RNAseq rules
    #"reports/16_gene_unity_setup.html",
    # Spatial rules
    # "src/spatial/organize_vectra_samples.html",
    # "reports/spatial/01_summary_QC.html",
    # expand("reports/spatial/object_qc_by_batch/02_object_QC_{batch}.html", batch = ['batch' + str(i) for i in range (1, 8)]),
    # expand("reports/spatial/batch_marker_viz/03_{batch}_marker_coexpression.html", batch = ['batch' + str(i) for i in range (1, 8)]),
    # "reports/spatial/03_aggregate_marker_combos.html",
    # "reports/spatial/04_report_cell_types.html",
    # "reports/spatial/05_density.html",
    # "reports/spatial/06_kruskal_total_density.html",
    # "reports/spatial/06b_inv_time_density.html",
    # expand("reports/spatial/07_cox_total_density_{outcome}.html", outcome = ['OS','DRS']),
    # "reports/spatial/08_ig_clusters_cd20.html",
    # "reports/spatial/09_tissue_segmentation.html"



