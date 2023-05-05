# Post-partum breast cancer

This repository contains the code for the analysis of RNA-seq data from the KWF project "Postpartum breast cancer diagnosed during involution: a distinct entity with unique clinicopathological, molecular and immunological features".

## Approach

This project utilizes a three-pronged approach:

1) Survival analysis based on clinical outcomes within the PPBC cohort
2) RNAseq-based analyses on FFPE-preserved primary tumors from a subset of patients within the PPBC cohort
3) Spatial analysis from Vectra multiplex panels on a subset of the aforementioned patients

## Instructions

1) Install R and Python and their dependencies. R dependency versions are in `renv.lock` and can be restored with `renv::restore()`. For Python, use `pipenv install` to generate a virtual environment based on `requirements.txt`, and `pipenv shell` to activate it.
2) Add the RNA fastq files to `data/rnaseq/RAW` and the spatial vectra files to `data/vectra/raw`.
3) Run with `snakemake`. If you get a permission error when knitting a notebook, run snakemake again. This is due to rmarkdown's inability to handle multiple simultaneous knitting jobs from the same file, and should resolve itself on subsequent runs.

Note: to get the figures numbered as they are in the published article, you will need to type `snakemake renumber_figs` into the console.
