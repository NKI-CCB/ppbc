library(here)
library(tidyverse)

appDir <- here::here("shinyApp", "VisualizePPBCgene")
stopifnot(dir.exists(appDir))
dataDir <-file.path(appDir, "data")
print(dataDir)
dir.create(dataDir, showWarnings = F)
stopifnot(dir.exists(dataDir))

#Gene annotation data
#file.copy(here("data","Rds","12e_gx_annot.Rds"), file.path(dataDir,"12e_gx_annot.Rds"), overwrite = T)
gx_annot <- readRDS(here::here("data","Rds","12e_gx_annot.Rds"))
more_ids <- read_tsv(here::here("data/external/ensembl_universal_ids_v94.txt"))

#Some uniprot ids are doubly presented with both text and NA, unclear why
head(arrange(more_ids, `Gene stable ID`))
nrow(more_ids) #38987

#Split up and remove NAs for each column of interest
uniprot_ids <- select(more_ids,
       ensembl_gene_id = `Gene stable ID`,
       uniprot_id = `UniProtKB/Swiss-Prot ID`) %>%
  filter(!is.na(uniprot_id)) %>%
  distinct() %>%
  arrange(ensembl_gene_id)
nrow(uniprot_ids) #21487
head(uniprot_ids)

entrez_ids <- select(more_ids,
                      ensembl_gene_id = `Gene stable ID`,
                      entrez_id = `NCBI gene ID`) %>%
  filter(!is.na(entrez_id)) %>%
  arrange(ensembl_gene_id) %>%
  distinct()
nrow(entrez_ids) #22570
head(entrez_ids)

gx_annot <- left_join(gx_annot,
                      uniprot_ids,
                      by = "ensembl_gene_id") %>%
  left_join(., entrez_ids, by = "ensembl_gene_id") %>%
  arrange(ensembl_gene_id)

gx_annot$gene_type <- str_replace_all(gx_annot$gene_type, "_", " ")
head(gx_annot)

saveRDS(gx_annot, file.path(dataDir,"12e_gx_annot.Rds"))

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

