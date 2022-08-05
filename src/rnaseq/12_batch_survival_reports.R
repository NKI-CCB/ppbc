library(here)
library(rmarkdown)
library(tidyverse)

rm(list = ls())

coxfiles = c(
  "12_multi_genewise_os.Rds", "12_uni_genewise_os.Rds",
  "12_multi_genewise_drs.Rds", "12_uni_genewise_drs.Rds",
  "12_inv_multi_genewise_os.Rds", "12_inv_uni_genewise_os.Rds",
  "12_inv_multi_genewise_drs.Rds", "12_inv_uni_genewise_drs.Rds"
)

projDir <- "/DATA/share/postpartumbc"
Rmd <- file.path(projDir, "reports/12_genewise_survival.Rmd")
stopifnot(file.exists(Rmd))
overwrite = F

#rmarkdown::render(p$rmd, rmarkdown::md_document(variant='markdown'), p$out, params=p$params)
for(cox in coxfiles){
  thisrds=file.path(projDir, "data", "Rds", cox)
  print(thisrds)
  stopifnot(file.exists(thisrds))
  html=file.path(projDir, "reports", str_replace(cox, "Rds", "html"))
  print(html)
  if (!file.exists(html) | overwrite == T){
    rmarkdown::render(input =  Rmd,
                      output_file = html,
                      params = list(survival_results = thisrds))
  }
                    
  
}


