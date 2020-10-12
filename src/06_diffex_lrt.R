#Open the bash terminal and log into harris/darwin
#Run this command:
#export OMP_NUM_THREADS=1
#Then run the script from the console

library(DESeq2)
library(apeglm)
library(here)
library(tidyverse)
library(tictoc)

#### Load and pre-process data ####

#sink("testingparalellel.txt")
rm(list = ls())
dds <- readRDS(file = here("data/Rds/05b_dds_filtered.Rds"))

print(design(dds))

gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>%
  select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>%
  distinct()

blind_vst = function(dds){
  design(dds) = formula(~ 1)
  vsd = vst(dds, blind=T)
  #mat = assay(vsd)
  return(vsd)
}

#Whether to overwrite pre-existing results files
overwrite <- T

#Create a version without the lactation group, as it has very few samples compared to the rest
dds.nolac = dds[,colnames(dds)[colData(dds)$study_group != "ppbc_lac"]]
dds.nolac$study_group = droplevels(dds.nolac$study_group)
dds.nolac$study_group = factor(dds.nolac$study_group,
                               c("non_prbc", "prbc", "ppbc_inv"))
dds.nolac$PPBC = droplevels(dds.nolac$PPBC)
dds.nolac$PPBC = factor(dds.nolac$PPBC,
                        c("nulliparous", "pregnant", "involuting"))

#### Variance stabilizing transformations ####

if(file.exists(here("data","Rds","06_vsd.Rds")) == F | overwrite == T){
  print("Variance stabilizing transformation")
  tic()
  vsd = blind_vst(dds)
  saveRDS(vsd, here("data","Rds","06_vsd.Rds"))
  toc()
}

if(file.exists(here("data","Rds","06_vsd_nolac.Rds")) == F | overwrite == T){
  print("Variance stabilizing transformation")
  tic()
  vsd.nolac = blind_vst(dds.nolac)
  #print(table(vsd.nolac$study_group))
  toc()
  saveRDS(vsd.nolac, here("data","Rds","06_vsd_nolac.Rds"))
}


#### Likelihood ratio tests ####

#The likelihood ratio test (LRT) is an ANOVA-like test design to check for genes which are significantly differentially expressed in at least one group. 
#From the manual: "The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model, in which some of the terms of the full model are removed.
#The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero."


#Our design formula is as follows:
print(design(dds))

### Study_group ####

#We will test it versus a formula that lacks the study group for genes of interest.
rds = here("data", "Rds", "06_ddsLRT.Rds")

if(file.exists(rds)==F | overwrite == T){
  print("Starting LRT vs reduced formula ~ batch + PAM50")
  tic()
  ddsLRT = DESeq(dds, test="LRT", reduced= ~ batch + PAM50)
  toc()
  saveRDS(object = ddsLRT, file = rds)
} 

#Again without lac

rds = here("data", "Rds", "06_ddsLRT_nolac.Rds")

if(file.exists(rds)==F | overwrite == T){
  print("Starting LRT vs reduced formula ~ batch + PAM50 without lac")
  tic()
  ddsLRT = DESeq(dds.nolac, test="LRT", reduced= ~ batch + PAM50)
  toc()
  saveRDS(object = ddsLRT, file = rds)
}

#### PAM50 ####

rds = here("data", "Rds", "06_ddsLRT_pam.Rds")

if(file.exists(rds)==F | overwrite == T){
  print("Starting LRT vs reduced formula ~ batch + study_group")
  tic()
  ddsLRT = DESeq(dds, test="LRT", reduced= ~ batch + study_group)
  toc()
  saveRDS(object = ddsLRT, file = rds)
}

#Again without lac

rds = here("data", "Rds", "06_ddsLRT_pam_nolac.Rds")

if(file.exists(rds)==F | overwrite == T){
  print("Starting LRT vs reduced formula ~ batch + study_group without lac")
  tic()
  ddsLRT = DESeq(dds.nolac, test="LRT", reduced= ~ batch + study_group)
  toc()
  saveRDS(object = ddsLRT, file = rds)
}

#### Batch ####

rds = here("data", "Rds", "06_ddsLRT_batch.Rds")

if(file.exists(rds)==F | overwrite == T){
  print("Starting LRT vs reduced formula ~ PAM50 + study_group")
  tic()
  ddsLRT = DESeq(dds.nolac, test="LRT", reduced= ~ PAM50 + study_group)
  toc()
  saveRDS(object = ddsLRT, file = rds)
}


#Again without lac

rds = here("data", "Rds", "06_ddsLRT_batch_nolac.Rds")

if(file.exists(rds)==F | overwrite == T){
  print("Starting LRT vs reduced formula ~ PAM50 + study_group without lac")
  tic()
  ddsLRT = DESeq(dds.nolac, test="LRT", reduced= ~ PAM50 + study_group)
  toc()
  saveRDS(object = ddsLRT, file = rds)
}
