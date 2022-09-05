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

rule all:
  input:
    # RNAseq rules
    expand("data/rnaseq/salmon/{sample}/quant.sf", sample=config['samples']),
    expand("results/rnaseq/fastqc/{sample}_fastqc.html", sample=config['samples']),
    "reports/rnaseq/08_diffex_onevsrest.html",
    "reports/rnaseq/09_diffex_time_involution.html",
    "reports/rnaseq/09b_diffex_time_breastfeeding.html",
    expand("reports/rnaseq/14_subgroup_diffex_{comp}.html",
      comp=["ppbcpw_vs_npbc","ppbcpw_vs_prbc","ppbcpw_vs_rest"]),
    "reports/rnaseq/10_CibersortX.pdf",
    "reports/rnaseq/11_clustering_involution.html",
    expand("reports/rnaseq/12_{cox}.html", cox=genewise_cox+interaction_cox),
    "results/rnaseq/TRUST/TRUST_results.xlsx",
    "reports/rnaseq/16_gene_unity_setup.html",
    # Spatial rules
    "reports/spatial/05_report_cell_types.html",
    "reports/spatial/07_kruskal_density.html",
    "reports/spatial/08_inv_time_density.html",
    expand("reports/spatial/09_cox_total_density_{outcome}.html", outcome = ['OS','DRS']),
    # "reports/spatial/10_ig_clusters_cd20.html", #Needs RNA input
    "reports/spatial/12_spatstat_overview.html"

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
