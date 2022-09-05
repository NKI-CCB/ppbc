library(here)
library(tidyverse)

appDir <- here::here("shinyApp", "VisualizePPBCgene")
stopifnot(dir.exists(appDir))
dataDir <-file.path(appDir, "data")
print(dataDir)
dir.create(dataDir, showWarnings = F)
stopifnot(dir.exists(dataDir))

#Need uniprot ids in app
#Gene annotation data
gx_annot <- readRDS(here("data/rnaseq/processed/bx_annot.Rds"))
saveRDS(gx_annot, file.path(dataDir,"app_gx_annot.Rds"))

#Genewise overall survival and drs
file.copy(here("results", "rnaseq", "survival", "12_cox_allgenes.xlsx"),
          file.path(dataDir, "12_cox_allgenes.xlsx"), overwrite = T)

#Differential expression results
lrt <- here("results", "rnaseq", "diffex", "06_LRT_allgenes.xlsx") %>% 
  readxl::excel_sheets() %>% 
  purrr::set_names() %>% 
  map(readxl::read_excel, path = here("results", "rnaseq", "diffex", "06_LRT_allgenes.xlsx"))
names(lrt)

pw <- here("results", "rnaseq", "diffex", "07_pairwise_comparisons_allgenes.xlsx") %>% 
  readxl::excel_sheets() %>% 
  purrr::set_names() %>% 
  map(readxl::read_excel, path = here("results", "rnaseq", "diffex", "07_pairwise_comparisons_allgenes.xlsx"))
names(pw)

ovr = here("results", "rnaseq", "diffex", "08_one_vs_rest_allgenes.xlsx") %>% 
  readxl::excel_sheets() %>% 
  purrr::set_names() %>% 
  map(readxl::read_excel, path = here("results", "rnaseq", "diffex", "08_one_vs_rest_allgenes.xlsx"))
names(ovr)

#Substitute lrt down to all genes and rename it for convenience.
lrt = lrt["repLRT_all"]
names(lrt) = "LRT"

#Combine all results into single list.
res_list = append(lrt, pw)
res_list = append(res_list, ovr)
names(res_list) = str_remove_all(str_remove_all(names(res_list), "rep_"), "_all")
names(res_list) = str_replace(names(res_list), "_", "_vs_")

saveRDS(res_list, file.path(dataDir, "app_diffex_res_list.Rds"))

#Sample metadata
coxdata = readRDS(here("data", "rnaseq", "processed", "12_coxdata.Rds"))
gene_col = which(colnames(coxdata)=="ENSG00000000003")
sample_data = coxdata[,1:(gene_col-1)]
#head(sample_data)
#Binary yes-no for involution
sample_data = sample_data %>%
  mutate(involution = if_else(study_group == "ppbcpw", 1, 0)) %>%
  select(sample_name:PPBC, involution, everything())

saveRDS(sample_data, file.path(dataDir,"app_survival_sample_data.Rds"))


#TMM/log normalized gene expression matrices, with ensembl ids and gene symbols
#ensembl ids:
ens_mat <- t(coxdata[gene_col:ncol(coxdata)])
#dim(ens_mat)
#ens_mat[1:10,1:10]
saveRDS(ens_mat, file.path(dataDir, "app_ensembl_tmmnorm_genesxsample.Rds"))

#gene symbols:
source(here("src/rnaseq/summarize_duplicate_ids.R"))
geneEx = rownames_to_column(as.data.frame(ens_mat), "ensembl_gene_id")
geneEx = right_join(select(gx_annot, gene_name, ensembl_gene_id),
                    geneEx, by = "ensembl_gene_id") %>%
  select(-ensembl_gene_id)
geneEx = summarize_expression_duplicate_ids(geneEx, id_column = "gene_name")
rownames(geneEx) = NULL
geneEx = column_to_rownames(geneEx, "GeneSymbol")
geneEx = as.matrix(geneEx)

saveRDS(geneEx, file.path(dataDir, "app_symbol_tmmnorm_genesxsample.Rds"))
