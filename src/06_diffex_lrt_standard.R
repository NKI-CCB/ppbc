#Open the bash terminal and log into harris/darwin
#Run this command:
#export OMP_NUM_THREADS=1
#Then run the script from the console

library("DESeq2")
library(apeglm)
library(here)
library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(openxlsx)
library(scrime)
library(RColorBrewer)

#### Load and pre-process data ####

#sink("testingparalellel.txt")

dds <- readRDS(file = here("data/Rds/05_dds_PAM50_batch.Rds"))

print(design(dds))

gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>% distinct()

source(here("src", "deseq_report_functions.R"))

#Whether to overwrite pre-existing results files
overwrite <- F

## Filtering

#Require non-zero counts in 1/3 of the dataset.
keep <- rowSums(counts(dds)!=0) >= ceiling(ncol(dds)/3)
table(keep)
dds <- dds[keep,]

#Create a version without the lactation group, as it has very few samples compared to the rest
dds.nolac = dds[,colnames(dds)[colData(dds)$study_group != "ppbc_lac"]]
dds.nolac$study_group = droplevels(dds.nolac$study_group)

#### Variance stabilizing transformations ####

if(file.exists(here("data","Rds","06_vsd_standardfilter.Rds")) == F & overwrite == T){
  print("Variance stabilizing transformation")
  vsd = blind_vst(dds)
  saveRDS(vsd, here("data","Rds","06_vsd_standardfilter.Rds"))
}

if(file.exists(here("data","Rds","06_vsd_nolac.Rds")) == F & overwrite == T){
  print("Variance stabilizing transformation")
  vsd.nolac = blind_vst(dds.nolac)
  print(table(vsd.nolac$study_group))
  saveRDS(vsd.nolac, here("data","Rds","06_vsd_nolac.Rds"))
}


#### Likelihood ratio tests ####

#The likelihood ratio test (LRT) is an ANOVA-like test design to check for genes which are significantly differentially expressed in at least one group. From the manual: "The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model, in which some of the terms of the full model are removed. The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero."
#Our design formula is as follows (tumor_purity has been removed on the basis of notebook 5):

print(design(dds))

### Study_group ####

#We will test it versus a formula that lacks the study group for genes of interest.
if(file.exists(here("data", "Rds", "06_ddsLRT_standard.Rds"))==F){
  print("Starting LRT vs reduced formula ~ batch + PAM50")
  start <- Sys.time()
  ddsLRT = DESeq(dds, test="LRT", reduced= ~ batch + PAM50)
  end <- Sys.time()
  print(end-start)
  saveRDS(ddsLRT, here("data", "Rds", "06_ddsLRT_standard.Rds"))
} else if(file.exists(here("data", "Rds", "06_ddsLRT_standard.Rds"))==T & overwrite == T){
  print("06_ddsLRT_standard.Rds exists and overwrite is true")
  print("Starting LRT vs reduced formula ~ batch + PAM50")
  start <- Sys.time()
  ddsLRT = DESeq(dds, test="LRT", reduced= ~ batch + PAM50)
  end <- Sys.time()
  print(end-start)
  saveRDS(ddsLRT, here("data", "Rds", "06_ddsLRT_standard.Rds"))
} else if(file.exists(here("data", "Rds", "06_ddsLRT_standard.Rds"))==T & overwrite == F){
  print("06_ddsLRT_standard.Rds exists and overwrite is false")
  print("Skipping to next")
}

#Again without lac
if(file.exists(here("data", "Rds", "06_ddsLRT_nolac.Rds"))==F){
  print("Starting LRT vs reduced formula ~ batch + PAM50 without lac")
  start <- Sys.time()
  ddsLRT.nolac = DESeq(dds.nolac, test="LRT", reduced= ~ batch + PAM50)
  end <- Sys.time()
  print(end-start)
  saveRDS(ddsLRT.nolac, here("data", "Rds", "06_ddsLRT_nolac.Rds"))
} else if(file.exists(here("data", "Rds", "06_ddsLRT_nolac.Rds"))==T & overwrite == T){
  print("06_ddsLRT_nolac.Rds exists and overwrite is true")
  print("Starting LRT vs reduced formula ~ batch + PAM50 without lac")
  start <- Sys.time()
  ddsLRT.nolac = DESeq(dds.nolac, test="LRT", reduced= ~ batch + PAM50)
  end <- Sys.time()
  print(end-start)
  saveRDS(ddsLRT.nolac,here("data", "Rds", "06_ddsLRT_nolac.Rds"))
} else if(file.exists(here("data", "Rds", "06_ddsLRT_nolac.Rds"))==T & overwrite == F){
  print("06_ddsLRT_nolac.Rds exists and overwrite is false")
  print("Skipping to next")
}

#### PAM50 ####

if(file.exists(here("data", "Rds", "06_ddsLRT_pam.Rds"))==F){
  print("Starting LRT vs reduced formula ~ batch + study_group")
  start <- Sys.time()
  ddsLRT = DESeq(dds, test="LRT", reduced= ~ batch + study_group)
  end <- Sys.time()
  print(end-start)
  saveRDS(ddsLRT, here("data", "Rds", "06_ddsLRT_pam.Rds"))
} else if(file.exists(here("data", "Rds", "06_ddsLRT_pam.Rds"))==T & overwrite == T){
  print("06_ddsLRT_pam.Rds exists and overwrite is true")
  print("Starting LRT vs reduced formula ~ batch + study_group")
  start <- Sys.time()
  ddsLRT = DESeq(dds, test="LRT", reduced= ~ batch + study_group)
  end <- Sys.time()
  print(end-start)
  saveRDS(ddsLRT, here("data", "Rds", "06_ddsLRT_pam.Rds"))
} else if(file.exists(here("data", "Rds", "06_ddsLRT_pam.Rds"))==T & overwrite == F){
  print("06_ddsLRT_pam.Rds exists and overwrite is false")
  print("Skipping to next")
}

#Again without lac
if(file.exists(here("data", "Rds", "06_ddsLRT_pam_nolac.Rds"))==F){
  print("Starting LRT vs reduced formula ~ batch + study_group without lac")
  start <- Sys.time()
  ddsLRT.nolac = DESeq(dds.nolac, test="LRT", reduced= ~ batch + study_group)
  end <- Sys.time()
  print(end-start)
  saveRDS(ddsLRT.nolac,here("data", "Rds", "06_ddsLRT_pam_nolac.Rds"))
} else if(file.exists(here("data", "Rds", "06_ddsLRT_pam_nolac.Rds"))==T & overwrite == T){
  print("06_ddsLRT_pam_nolac.Rds exists and overwrite is true")
  print("Starting LRT vs reduced formula ~ batch + study_group without lac")
  start <- Sys.time()
  ddsLRT.nolac = DESeq(dds.nolac, test="LRT", reduced= ~ batch + study_group)
  end <- Sys.time()
  print(end-start)
  saveRDS(ddsLRT.nolac,here("data", "Rds", "06_ddsLRT_pam_nolac.Rds"))
} else if(file.exists(here("data", "Rds", "06_ddsLRT_pam_nolac.Rds"))==T & overwrite == F){
  print("06_ddsLRT_pam_nolac.Rds exists and overwrite is false")
  print("Skipping to next")
}

#### Batch ####

if(file.exists(here("data", "Rds", "06_ddsLRT_batch.Rds"))==F){
  print("Starting LRT vs reduced formula ~ PAM50 + study_group")
  start <- Sys.time()
  ddsLRT = DESeq(dds, test="LRT", reduced= ~ PAM50 + study_group)
  end <- Sys.time()
  print(end-start)
  saveRDS(ddsLRT, here("data", "Rds", "06_ddsLRT_batch.Rds"))
} else if(file.exists(here("data", "Rds", "06_ddsLRT_batch.Rds"))==T & overwrite == T){
  print("06_ddsLRT_batch.Rds exists and overwrite is true")
  print("Starting LRT vs reduced formula ~ PAM50 + study_group")
  start <- Sys.time()
  ddsLRT = DESeq(dds, test="LRT", reduced= ~ PAM50 + study_group)
  end <- Sys.time()
  print(end-start)
  saveRDS(ddsLRT, here("data", "Rds", "06_ddsLRT_batch.Rds"))
} else if(file.exists(here("data", "Rds", "06_ddsLRT_batch.Rds"))==T & overwrite == F){
  print("06_ddsLRT_batch.Rds exists and overwrite is false")
  print("Skipping to next")
}

#Again without lac
if(file.exists(here("data", "Rds", "06_ddsLRT_batch_nolac.Rds"))==F){
  print("Starting LRT vs reduced formula ~ PAM50 + study_group without lac")
  start <- Sys.time()
  ddsLRT.nolac = DESeq(dds.nolac, test="LRT", reduced= ~ PAM50 + study_group)
  end <- Sys.time()
  print(end-start)
  saveRDS(ddsLRT.nolac,here("data", "Rds", "06_ddsLRT_batch_nolac.Rds"))
} else if(file.exists(here("data", "Rds", "06_ddsLRT_batch_nolac.Rds"))==T & overwrite == T){
  print("06_ddsLRT_batch_nolac.Rds exists and overwrite is true")
  print("Starting LRT vs reduced formula ~ PAM50 + study_group without lac")
  start <- Sys.time()
  ddsLRT.nolac = DESeq(dds.nolac, test="LRT", reduced= ~ PAM50 + study_group)
  end <- Sys.time()
  print(end-start)
  saveRDS(ddsLRT.nolac,here("data", "Rds", "06_ddsLRT_batch_nolac.Rds"))
} else if(file.exists(here("data", "Rds", "06_ddsLRT_batch_nolac.Rds"))==T & overwrite == F){
  print("06_ddsLRT_batch_nolac.Rds exists and overwrite is false")
  print("Skipping to next")
}