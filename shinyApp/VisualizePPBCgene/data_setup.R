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
#file.copy(here("data","Rds","12e_gx_annot.Rds"), file.path(dataDir,"12e_gx_annot.Rds"), overwrite = T)
gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>%
  select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>%
  distinct() %>%
  mutate(description = str_remove_all(description, " \\[.*\\]"))

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

saveRDS(gx_annot, file.path(dataDir,"app_gx_annot.Rds"))

#Genewise overall survival and drs
file.copy(here("results", "survival", "12_cox_allgenes.xlsx"), file.path(dataDir, "12_cox_allgenes.xlsx"), overwrite = T)

#Interaction model overall survival
file.copy(here("results", "survival", "13_multi_interaction_os.csv"), file.path(dataDir, "13_multi_interaction_os.csv"), overwrite = T)

#Interaction model distant recurrence
file.copy(here("results", "survival", "13_multi_interaction_drs.csv"), file.path(dataDir, "13_multi_interaction_drs.csv"), overwrite = T)

#Differential expression results
lrt <- here("results", "diffex", "06_LRT_allgenes.xlsx") %>% 
  readxl::excel_sheets() %>% 
  purrr::set_names() %>% 
  map(readxl::read_excel, path = here("results", "diffex", "06_LRT_allgenes.xlsx"))
names(lrt)

pw <- here("results", "diffex", "07_pairwise_comparisons_allgenes.xlsx") %>% 
  readxl::excel_sheets() %>% 
  purrr::set_names() %>% 
  map(readxl::read_excel, path = here("results", "diffex", "07_pairwise_comparisons_allgenes.xlsx"))
names(pw)

ovr = here("results", "diffex", "08_one_vs_rest_allgenes.xlsx") %>% 
  readxl::excel_sheets() %>% 
  purrr::set_names() %>% 
  map(readxl::read_excel, path = here("results", "diffex", "08_one_vs_rest_allgenes.xlsx"))
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
coxdata = readRDS(here("data", "Rds", "12_coxdata.Rds"))
gene_col = which(colnames(coxdata)=="ENSG00000000003")
sample_data = coxdata[,1:(gene_col-1)]
#head(sample_data)
#Binary yes-no for involution
sample_data = sample_data %>%
  mutate(involution = if_else(study_group == "ppbc_inv", 1, 0)) %>%
  select(sample_name:PPBC, involution, everything())

saveRDS(sample_data, file.path(dataDir,"app_survival_sample_data.Rds"))


#TMM/log normalized gene expression matrices, with ensembl ids and gene symbols
#ensembl ids:
ens_mat <- t(coxdata[gene_col:ncol(coxdata)])
#dim(ens_mat)
#ens_mat[1:10,1:10]
saveRDS(ens_mat, file.path(dataDir, "app_ensembl_tmmnorm_genesxsample.Rds"))

#gene symbols:
source(here("src", "general_R_tools.R"))
geneEx = rownames_to_column(as.data.frame(ens_mat), "ensembl_gene_id")
geneEx = right_join(select(gx_annot, gene_name, ensembl_gene_id),
                    geneEx, by = "ensembl_gene_id") %>%
  select(-ensembl_gene_id)
geneEx = summarize_expression_duplicate_ids(geneEx, id_column = "gene_name")
rownames(geneEx) = NULL
geneEx = column_to_rownames(geneEx, "GeneSymbol")
geneEx = as.matrix(geneEx)

saveRDS(geneEx, file.path(dataDir, "app_symbol_tmmnorm_genesxsample.Rds"))


