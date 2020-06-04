library(here)
library(rmarkdown)
library(tidyverse)

rm(list = ls())

genes_to_report = read_csv(
  here("reports", "genes_to_report")
)

projDir <- "/DATA/share/postpartumbc"
Rmd <- file.path(projDir, "reports/17_gene_report_template.Rmd")
stopifnot(file.exists(Rmd))
overwrite = F

outdir=here("reports/gene_reports")
dir.create(outdir, showWarnings = F)

for(gene in genes_to_report$Gene){
  print(gene)
  
  outfile=file.path(outdir, paste0(gene, "_report.html"))
  
  if (!file.exists(outfile) | overwrite == T){
    rmarkdown::render(input =  Rmd,
                      output_file = outfile,
                      params = list(gene = gene,
                                    show_functions = F,
                                    show_code = T))
  }
}


