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
    expand("data/rnaseq/salmon/{sample}/quant.sf", sample=config['samples']),
    #"reports/16_gene_unity_setup.html",
    # Spatial rules
    "reports/spatial/05_report_cell_types.html",
    "reports/spatial/06_density.html",
    "reports/spatial/07_kruskal_density.html",
    # "reports/spatial/08_inv_time_density.html",
    # expand("reports/spatial/09_cox_total_density_{outcome}.html", outcome = ['OS','DRS']),
    # "reports/spatial/10_ig_clusters_cd20.html", #Needs RNA input
    # "reports/spatial/12_tissue_segmentation.html"

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
