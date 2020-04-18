requireNamespace('rmarkdown')

p <- local({
  p <- list()
  args <- commandArgs(T)
  stopifnot(length(args) >= 2)
  p$rmd <- args[1]
  p$out  <- args[2]
  
  if (length(args) > 2) {
    params <- args[3:length(args)]
    stopifnot(length(params) %% 2 == 0)
    params_names <- params[seq(1, length(params), 2)]
    params_values <- params[seq(1, length(params), 2)+1]
    stopifnot(stringr::str_sub(params_names, 1, 2) == '--')
    p$params <- structure(as.list(params_values),
                          names=stringr::str_sub(params_names, 3, -1))
  } else {
    p$params=NULL
  }
  p
})

#rmarkdown::render(p$rmd, rmarkdown::md_document(variant='markdown'), p$out, params=p$params)
rmarkdown::render(input = p$rmd, output_file = p$out, params=p$params)

#Example usage of the rmarkdown script with Snakemake
#rule knit_pathway_analysis_mouse:
#  output:
#  html="reports/pathway-analysis-nki_fsh2-{treatment}.html",
#input:
#  rmd="src/reports/pathway-analysisi-mouse.Rmd",
#  script="src/reports/rmarkdown.R",
#  rnaseq="data/processed/gene-expression-nki_fsh2.tsv",
#  mgi_id_map="data/external/gencode/ensembl_mgi.tsv",
#  mgi_homologs="data/external/mgi/HOM_AllOrganism.rpt",
#  sample_sheet="data/raw/sample_annotation/sample_sheet_nki_fsh2.tsv",
#  gene_sets_c2_cgp="data/external/msigdb/c2.cgp.v6.2.entrez.gmt",
#  gene_sets_c2_cp="data/external/msigdb/c2.cp.v6.2.entrez.gmt",
#  gene_sets_h_all="data/external/msigdb/h.all.v6.2.entrez.gmt",
#wildcard_constraints:
#  treatment='FSHu|FSHg'
#shell:
#  "mkdir -p reports\n"
#  "Rscript {input.script} {input.rmd} $PWD/{output.html} "
#  "--rnaseq {input.rnaseq} "
#  "--mgi_id_map {input.mgi_id_map} "
#  "--mgi_homologs {input.mgi_homologs} "
#  "--sample_sheet {input.sample_sheet} "
#  "--gene_sets data/external/msigdb "
#  "--treatment {wildcards.treatment}"


