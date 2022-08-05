library(here)
library(rmarkdown)
library(tidyverse)

rm(list = ls())

gene = "UBD"

projDir <- "/DATA/share/postpartumbc"
Rmd <- file.path(projDir, "reports/17_gene_report_template.Rmd")
stopifnot(file.exists(Rmd))
overwrite = T

outdir=here("reports/gene_reports")
dir.create(outdir, showWarnings = F)
outfile=file.path(outdir, paste0(gene, "_report.html"))

if (!file.exists(outfile) | overwrite == T){
  rmarkdown::render(input =  Rmd,
                    output_file = outfile,
                    params = list(gene = gene,
                                  show_functions = F,
                                  show_code = T))
}

