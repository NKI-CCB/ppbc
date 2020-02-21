library(here)
library(tidyverse)

appDir <- here("shinyApp", "VisualizePPBCgene")
stopifnot(dir.exists(appDir))
dataDir <-file.path(appDir, "data")
print(dataDir)
dir.create(dataDir, showWarnings = F)
stopifnot(dir.exists(dataDir))

#Gene annotation data
file.copy(here("data","Rds","12e_gx_annot.Rds"), file.path(dataDir,"12e_gx_annot.Rds"), overwrite = T)

#Genewise overall survival and drs
file.copy(here("results", "survival", "12_cox_allgenes.xlsx"), file.path(dataDir, "12_cox_allgenes.xlsx"), overwrite = T)

#Interaction model overall survival
file.copy(here("results", "survival", "12d_os_gene:inv_newformula.xlsx"), file.path(dataDir, "12d_os_gene:inv_newformula.xlsx"), overwrite = T)

#Interaction model distant recurrence
file.copy(here("results", "survival", "12d_drs_gene:inv_newformula.xlsx"), file.path(dataDir, "12d_drs_gene:inv_newformula.xlsx"), overwrite = T)

#Differential expression results
file.copy(here("data", "Rds", "12e_diffex_res_list.Rds"), file.path(dataDir, "12e_diffex_res_list.Rds"), overwrite = T)

#Sample metadata
file.copy(here("data","Rds","12e_survival_sample_data.Rds"), file.path(dataDir, "12e_survival_sample_data.Rds"), overwrite = T)

#TMM/log normalized gene expression matrices, with ensembl ids and gene symbols
file.copy(here("data", "Rds", "12e_symbol_tmmnorm_genesxsample.Rds"), file.path(dataDir, "12e_symbol_tmmnorm_genesxsample.Rds"), overwrite = T)

file.copy(here("data", "Rds", "12e_ensembl_tmmnorm_genesxsample.Rds"), file.path(dataDir, "12e_ensembl_tmmnorm_genesxsample.Rds"), overwrite = T)

