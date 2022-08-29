library(survival)
library(survminer)
library(RColorBrewer)
library(glmnet)
library(tidyverse)
library(here)

overwrite <- T
set.seed(123)

resDir = here("data/rnaseq/processed")
dir.create(resDir, showWarnings = F)
stopifnot(file.exists(resDir))
source("src/rnaseq/genewise_survival.R")

#### Load data ----

invdata <- readRDS(here("data/rnaseq/processed/12_invdata.Rds"))
gx_annot <- read_tsv(here("data/rnaseq/metadata/01_gene_annot.tsv")) %>% 
  select(ensembl_gene_id = gene_id,
         gene_name, gene_type,
         description = gene_description) %>% distinct()

gene_col = which(colnames(invdata)=="ENSG00000000003")
print(paste("Gene columns from", gene_col, "to", ncol(invdata), "of invdata: involution only"))
print(head(invdata[,(gene_col-3):(gene_col+1)]))

#### Overall survival ----

#Testing
#testdata = invdata[,1:30]
#genewise_cox(
#  gene_list = colnames(testdata)[gene_col:ncol(testdata)],
#  data = testdata,
#  show_runtime = T,
#  time="time_OS_months",
#  event="death",
#  type="multivariate")

resPath = file.path(resDir, "12_inv_multi_genewise_os.Rds")

if(file.exists(resPath) == F| overwrite == T){
  print("Multivariate survival, involution samples")
  res <- genewise_cox(
    gene_list = colnames(invdata)[gene_col:ncol(invdata)],
    data = invdata,
    show_runtime = T,
    time="time_OS_months",
    event="death",
    type="multivariate")
  saveRDS(res, resPath)
}

#Univariate
resPath = file.path(resDir, "12_inv_uni_genewise_os.Rds")

#Testing
#testdata = invdata[,1:30]
#genewise_cox(
#  gene_list = colnames(testdata)[gene_col:ncol(testdata)],
#  data = testdata,
#  show_runtime = T,
#  time="time_OS_months",
#  event="death",
#  type="univariate")

if(file.exists(resPath) == F| overwrite == T){
  print("Univariate survival, involution samples")
  res <- genewise_cox(
    gene_list = colnames(invdata)[gene_col:ncol(invdata)],
    data = invdata,
    show_runtime = T,
    time="time_OS_months",
    event="death",
    type="univariate")
  saveRDS(res, resPath)
}

#### DRS ----

#Multivariate
resPath = file.path(resDir, "12_inv_multi_genewise_drs.Rds")

#genewise_cox(
#  gene_list = colnames(testdata)[gene_col:ncol(testdata)],
#  data = testdata,
#  show_runtime = T,
#  time="time_DRS_months",
#  event="distant_recurrence",
#  type="multivariate")

if(file.exists(resPath) == F| overwrite == T){
  print("Multivariate metastasis, involution samples")
  res <- genewise_cox(
    gene_list = colnames(invdata)[gene_col:ncol(invdata)],
    data = invdata,
    show_runtime = T,
    time="time_DRS_months",
    event="distant_recurrence",
    type="multivariate")
  saveRDS(res, resPath)
}

#Univariate
resPath = file.path(resDir, "12_inv_uni_genewise_drs.Rds")

#genewise_cox(
#  gene_list = colnames(testdata)[gene_col:ncol(testdata)],
#  data = testdata,
#  show_runtime = T,
#  time="time_DRS_months",
#  event="distant_recurrence",
#  type="univariate")

if(file.exists(resPath) == F| overwrite == T){
  print("Univariate metastasis, involution samples")
  res <- genewise_cox(
    gene_list = colnames(invdata)[gene_col:ncol(invdata)],
    data = invdata,
    show_runtime = T,
    time="time_DRS_months",
    event="distant_recurrence",
    type="univariate")
  saveRDS(res, resPath)
}
