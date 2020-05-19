library(here)
library(rmarkdown)
library(tidyverse)

rm(list = ls())

comps = c(
  "ppbc_inv_vs_non_prbc",
  "ppbc_inv_vs_prbc",
  "ppbc_inv_vs_rest"
)

projDir <- "/DATA/share/postpartumbc"
Rmd <- file.path(projDir, "reports/14_subgroup_diffex_by_comparison.Rmd")
stopifnot(file.exists(Rmd))
overwrite = F

#rmarkdown::render(p$rmd, rmarkdown::md_document(variant='markdown'), p$out, params=p$params)
for(comp in comps){
  html=paste0("14_subgroup_diffex_", comp, ".html")
  print(html)
  if (!file.exists(html) | overwrite == T){
    rmarkdown::render(input =  Rmd,
                      output_file = html,
                      params = list(comparison = comp))
  }
                    
  
}


