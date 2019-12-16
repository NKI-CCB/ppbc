library(DESeq2)
library(apeglm)
library(here)
library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(openxlsx)
library(scrime)
library(RColorBrewer)

#For staining, we have better access to Belgian patients than any other location.
#Check to see whether the IG signature is still present in the Belgium-only samples.

overwrite <- FALSE

#To speed up DESeq
#Open the bash terminal and log into harris/darwin
#Run this command:
#export OMP_NUM_THREADS=1
#Then run the script from the console

#### Load and pre-process data ####
dds <- readRDS(file = here("data/Rds/05_dds_PAM50_batch.Rds"))

gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>% distinct()

#### Filtering ####

#Minimum threshold is a nonzero count in at least a 3rd of all samples
keep <- rowSums(counts(dds)!=0) >= ceiling(ncol(dds)/3)
table(keep)

dds <- dds[keep,]

#Remove non Belgian samples
be_samples = as.data.frame(colData(dds)) %>%
  filter(country == "BE") %>%
  pull(sample_name)

print(paste("All samples: ", ncol(dds)))
dds <- dds[,be_samples]
print(paste("Belgian samples: ", ncol(dds)))

stopifnot(all(dds$country=="BE"))

table(dds$study_group)

## Retrieve immune genes

#For downstream analyses, we'd like to know which of the differentially expressed genes are immunologically relevant.
#An immune gene is defined as follows:

#Either immune/immuno/interleukin is part of the gene name and description OR the gene is part of the the ImmPort database.
#See https://www.innatedb.com/redirect.do?go=resourcesGeneLists

#The Immunology Database and Analysis Portal (ImmPort) system was developed under the Bioinformatics Integration Support Contract (BISC) Phase II by the Northrop Grumman Information Technology Health Solutions team for the NIH, NIAID, and DAIT. The principal investigator of the BISC project is Dr. Richard Scheuermann at University of Texas Southwestern Medical Center. The list of immunologically related genes in ImmPort is a collection of ~6,000 human genes, which was formed with the goal of retrieving all genes that have immune system-related functions. This list was generated using automatic searches of EntrezGene and Gene Ontology records using immunology-related keywords. The list was then manually curated by immunology experts examining various literature sources.

immune_gene_list <- read_csv(here("data", "external","gene-sets","InnateDB_genes.csv"))
head(immune_gene_list)

## Helper functions


source(here("src", "deseq_report_functions.R"))



#### One vs rest formula levels ####
print(design(dds))

dds$inv_vs_rest  = colData(dds) %>% as.data.frame() %>% mutate(inv_vs_rest = if_else(study_group == "ppbc_inv", "ppbc_inv", "rest")) %>%
  pull(inv_vs_rest) %>% factor(levels=c("rest", "ppbc_inv"))

dds$prbc_vs_rest  = colData(dds) %>% as.data.frame() %>% mutate(prbc_vs_rest = if_else(study_group == "prbc", "prbc", "rest")) %>%
  pull(prbc_vs_rest) %>% factor(levels=c("rest", "prbc"))

dds$lac_vs_rest  = colData(dds) %>% as.data.frame() %>% mutate(lac_vs_rest = if_else(study_group == "ppbc_lac", "ppbc_lac", "rest")) %>%
  pull(lac_vs_rest) %>% factor(levels=c("rest", "ppbc_lac"))

dds$nonprbc_vs_rest  = colData(dds) %>% as.data.frame() %>% mutate(nonprbc_vs_rest = if_else(study_group == "non_prbc", "non_prbc", "rest")) %>%
  pull(nonprbc_vs_rest) %>% factor(levels=c("rest", "non_prbc"))

table(dds$inv_vs_rest, dds$study_group) %>% print()
table(dds$prbc_vs_rest, dds$study_group) %>% print()
table(dds$lac_vs_rest, dds$study_group) %>% print()
table(dds$nonprbc_vs_rest, dds$study_group) %>% print()



#### Diffex: involution vs rest ####

dds$group = dds$inv_vs_rest

print(levels(dds$group))

design(dds) = ~batch + PAM50 + group

print(design(dds))

#Perform differential expression analysis.
resFile = here("data/Rds/13_dds_inv_vs_rest_BE.Rds")
if(file.exists(resFile)==F | overwrite == T){
  start = Sys.time()
  dds <- DESeq(dds)
  end = Sys.time()
  print(end-start)
  saveRDS(dds,resFile)
}

dds=readRDS(resFile)
levels(dds$group) %>% print()

#Shrink fold changes
resFile = here("data","Rds", "13_ape_inv_rest_BE.Rds")
if(file.exists(resFile)==F | overwrite == T){
  start = Sys.time()
  ape_inv_rest = shrinkRes(dds,  contrast=c("group", "ppbc_inv", "rest"))
  saveRDS(ape_inv_rest, resFile)
  end = Sys.time()
  print(end-start)
}