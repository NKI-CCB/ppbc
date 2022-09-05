library(here)
library(tidyverse)

gx_annot <- read_tsv(here("data/rnaseq/metadata/01_gene_annot.tsv"), show_col_types = F)
gx_annot = gx_annot %>%
  select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>%
  distinct() %>%
  mutate(description = str_remove_all(description, " \\[.*\\]"))

more_ids <- read_tsv(here::here("data/external/gene_ref/ensembl_universal_ids_v94.txt"),
                     show_col_types = F)

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

saveRDS(gx_annot, here("data/rnaseq/processed/bx_annot.Rds"))
