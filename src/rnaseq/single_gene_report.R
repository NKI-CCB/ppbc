library(here)
library(rmarkdown)
library(tidyverse)

gene = "UBD"
print(paste("Generating gene report for", gene))

projDir <- here()
Rmd <- file.path(projDir, "reports/rnaseq/17_gene_report_template.Rmd")
stopifnot(file.exists(Rmd))
overwrite = T

outdir=here("reports/rnaseq/gene_reports")
dir.create(outdir, showWarnings = F)
outfile=file.path(outdir, paste0(gene, "_report.html"))
print(paste("Outfile is:", outfile))

if (!file.exists(outfile) | overwrite == T){
  rmarkdown::render(input =  Rmd,
                    output_file = outfile,
                    params = list(gene = gene,
                                  show_functions = F,
                                  show_code = T))
}
