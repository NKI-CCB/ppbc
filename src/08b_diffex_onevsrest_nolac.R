library(DESeq2)
library(apeglm)
library(here)
library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(openxlsx)
library(scrime)
library(RColorBrewer)

#As notebook 8b
#To speed up DESeq
#Open the bash terminal and log into harris/darwin
#Run this command:
#export OMP_NUM_THREADS=1
#Then run the script from the console

overwrite <- FALSE
source(here("src", "deseq_report_functions.R"))

# Load and pre-process data
dds.nolac <- readRDS(file = here("data/Rds/06_dds_nolac.Rds"))

print("dds is pre-filtered")
print(nrow(dds.nolac))

# One vs rest formula levels
print(design(dds.nolac))

dds.nolac$inv_vs_rest  = colData(dds.nolac) %>% as.data.frame() %>% mutate(inv_vs_rest = if_else(study_group == "ppbc_inv", "ppbc_inv", "rest")) %>%
  pull(inv_vs_rest) %>% factor(levels=c("rest", "ppbc_inv"))

dds.nolac$prbc_vs_rest  = colData(dds.nolac) %>% as.data.frame() %>% mutate(prbc_vs_rest = if_else(study_group == "prbc", "prbc", "rest")) %>%
  pull(prbc_vs_rest) %>% factor(levels=c("rest", "prbc"))

dds.nolac$nonprbc_vs_rest  = colData(dds.nolac) %>% as.data.frame() %>% mutate(nonprbc_vs_rest = if_else(study_group == "non_prbc", "non_prbc", "rest")) %>%
  pull(nonprbc_vs_rest) %>% factor(levels=c("rest", "non_prbc"))

table(dds.nolac$inv_vs_rest, dds.nolac$study_group) %>% print()
table(dds.nolac$prbc_vs_rest, dds.nolac$study_group) %>% print()
table(dds.nolac$nonprbc_vs_rest, dds.nolac$study_group) %>% print()



#### Diffex: involution vs rest ####

#In this call to DESeq, we're focusing one one vs rest contrasts and using a Wald test instead of a likelihood ratio test.
print("Inv vs rest")
dds.nolac$group = dds.nolac$inv_vs_rest

print(levels(dds.nolac$group))

design(dds.nolac) = ~batch + PAM50 + group

print(design(dds.nolac))

#Perform differential expression analysis.
resFile = here("data/Rds/08b_dds_inv_vs_rest_nolac.Rds")
if(file.exists(resFile)==F | overwrite == T){
  start = Sys.time()
  dds.nolac <- DESeq(dds.nolac)
  end = Sys.time()
  print(end-start)
  saveRDS(dds.nolac,resFile)
}

dds.nolac=readRDS(resFile)
levels(dds.nolac$group) %>% print()

#Shrink fold changes
resFile = here("data","Rds", "08b_ape_inv_rest_nolac.Rds")
if(file.exists(resFile)==F | overwrite == T){
  start = Sys.time()
  ape_inv_rest_nolac = shrinkRes(dds.nolac,  contrast=c("group", "ppbc_inv", "rest"))
  saveRDS(ape_inv_rest_nolac, resFile)
  end = Sys.time()
  print(end-start)
}

#### Diffex: Pregnancy vs rest ####
dds.nolac = readRDS("data/Rds/08b_dds_inv_vs_rest_nolac.Rds")

print("Pregnancy vs rest")
dds.nolac$group = dds.nolac$prbc_vs_rest
levels(dds.nolac$group) %>% print()

#Perform differential expression analysis.

#From the manual:
#Although apeglm cannot be used with contrast, we note that many designs can be easily rearranged such that what was a contrast becomes its own coefficient.
#In this case, the dispersion does not have to be estimated again, as the designs are equivalent, up to the meaning of the coefficients.
#Instead, one need only run nbinomWaldTest to re-estimate MLE coefficients – these are necessary for apeglm – and then run lfcShrink specifying the coefficient of
#interest in resultsNames(dds.nolac). We give some examples below of producing equivalent designs for use with coef.
#We show how the coefficients change with model.matrix, but the user would, for example, either change the levels of dds$condition or replace the design using design(dds)<-,
#then run nbinomWaldTest followed by lfcShrink.



resFile = here("data/Rds/08b_dds_prbc_vs_rest_nolac.Rds")
if(file.exists(resFile)==F | overwrite == T){
  start = Sys.time()
  dds.nolac <- nbinomWaldTest(dds.nolac)
  end = Sys.time()
  print(end-start)
  saveRDS(dds.nolac, resFile)
}

dds.nolac=readRDS(resFile)
levels(dds.nolac$group) %>% print()

#Shrink fold changes
resFile = here("data","Rds", "08b_ape_prbc_rest_nolac.Rds")
if(file.exists(resFile)==F | overwrite == T){
  start = Sys.time()
  ape_prbc_rest_nolac = shrinkRes(dds.nolac,  contrast=c("group", "prbc", "rest"))
  saveRDS(ape_prbc_rest_nolac, resFile)
  end = Sys.time()
  print(end-start)
}

#### Diffex: Non-prbc vs rest ####
print("Nonprbc vs rest")
dds.nolac = readRDS("data/Rds/08b_dds_inv_vs_rest_nolac.Rds")

dds.nolac$group = dds.nolac$nonprbc_vs_rest
print(levels(dds.nolac$group))


#Perform differential expression analysis.
resFile = here("data/Rds/08b_dds_nonprbc_vs_rest_nolac.Rds")
if(file.exists(resFile)==F | overwrite == T){
  start = Sys.time()
  dds.nolac <- nbinomWaldTest(dds.nolac)
  end = Sys.time()
  print(end-start)
  saveRDS(dds.nolac, resFile)
}

dds.nolac=readRDS(resFile)
levels(dds.nolac$group) %>% print()

#Shrink fold changes
resFile = here("data","Rds", "08b_ape_nonprbc_rest_nolac.Rds")
if(file.exists(resFile)==F | overwrite == T){
  start = Sys.time()
  ape_nonprbc_rest_nolac = shrinkRes(dds.nolac,  contrast=c("group", "non_prbc", "rest"))
  saveRDS(ape_nonprbc_rest_nolac, resFile)
  end = Sys.time()
  print(end-start)
}
