rm(list = ls())

#devtools::install_github('kevinblighe/RegParallel')
#source("https://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#library(GEOquery) #Example only
#library(Biobase) #Example only
#library(RegParallel) #Example only
library(DESeq2)
library(survival)
library(survminer)
library(pheatmap)
library(RColorBrewer)
library(glmnet)
library(tidyverse)
library(here)
library(ggpubr)

#Whether to overwrite existing results files

overwrite <- FALSE

set.seed(123)

#### Set up data ####

# Pre-process data data
dds = readRDS(here("data/Rds/08_dds_inv_vs_rest_standard.Rds"))
gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>% distinct()


#Create data matrix with covariates from the sample data in the first few columns, then genes. Samples are rows.

## Prepare metadata

#Covariates: grade, stage, molecular subtype, age at diagnosis and year of diagnosis (Also possibly therapy type)
sd = as.data.frame(colData(dds)) %>%
  select(sample_name, death,
         years_overall_survival,
         months_of_followup,
         study_group, PAM50,
         grade, stage,
         inv_vs_rest, prbc_vs_rest,
         lac_vs_rest, nonprbc_vs_rest,
         year_of_diagnosis, age_at_diagnosis,
         overall_survival, months_overall_survival,
         distant_recurrence, months_to_drs,
         local_recurrence_free_survival, months_to_lrs)

#### Filter ####

#Exclude samples for which there has been less than 20 months of follow up, unless the reason for the lack of follow up is that the patient passed away due to disease (and not an accident or suicide).

#follow_up > 20 | death == 1

print(paste("Started with", nrow(sd), "samples"))

print("Samples lost due to lack of survival information:")
nrow(sd %>% filter(is.na(death)))
excluded = sd %>% filter(is.na(death)) %>%
  mutate(reason = "no_survival_data")

print("Samples excluded due to patient death unrelated to disease")
print(nrow(sd %>% filter(! death %in% c(NA, "due to disease", "no"))))
excluded = bind_rows(excluded,
                     mutate(filter(sd, !death %in% c(NA, "due to disease", "no")),
                            reason = "unrelated_death")
)
print("Samples due to insufficient follow up")
nrow(sd %>% filter(months_of_followup < 20 & death != "due to disease"))

excluded = bind_rows(excluded,
                     mutate(filter(sd, months_of_followup < 20 & !death == "due to disease"),
                            reason = "insufficient_followup"))

print("Samples with 999 (unknown) years of overall survival")
nrow(sd %>% filter(years_overall_survival == 999))

excluded = bind_rows(excluded,
                     mutate(filter(sd, years_overall_survival == 999),
                            reason = "unknown_survival_years"))

metadata = sd %>% filter(!is.na(death)) %>%
  filter(death != "other cause") %>%
  filter(years_overall_survival != 999) %>%
  filter(months_of_followup >= 20 | death == "due to disease")

print(paste(nrow(metadata), "samples remaining post filtering"))

stopifnot(nrow(filter(metadata,rowSums(is.na(metadata)) > 0)) == 0)

#### Transform gene data ####

## Prepare gene data
#fpm (fragments per million) is equivalent to cpm from edgeR
mat = t(log2(fpm(dds, robust=T) + 0.5))

#Exclude those samples which were also excluded from the metadata prep
mat = mat[rownames(mat) %in% metadata$sample_name,]
stopifnot(identical(rownames(mat), metadata$sample_name))

coxdata = cbind(metadata, mat)
#head(coxdata)

met_at_diag = colData(dds) %>% as.data.frame() %>%
  filter(months_to_drs == 0) %>% select(sample_name, months_to_drs, metastasis_at_diagnosis) %>%
  pull(sample_name)

#### Cox regression ####

#Results directory
resDir = here("data", "Rds")
dir.create(resDir, showWarnings = F)
stopifnot(file.exists(resDir))

#### Overall survival ####

#Multivariate
if(file.exists(file.path(resDir, "12_allgenes_multivar_surv.Rds")) == F | overwrite == T){
  print("Starting multivariate analysis for overall survival")
  start <- Sys.time()
  multivar_gene_formulas <- sapply(
    colnames(coxdata)[21:ncol(coxdata)],
    function(x) as.formula(
      paste(
        'Surv(time=months_overall_survival, event=overall_survival)~age_at_diagnosis+year_of_diagnosis+stage+grade+PAM50+',
        x)))

  multivar_gene_models <- lapply(multivar_gene_formulas, function(x){coxph(x, data = coxdata)})
  end <- Sys.time()
  print(end-start)

  saveRDS(multivar_gene_models, file=file.path(resDir, "12_allgenes_multivar_surv.Rds"))
}

#Univariate
if(file.exists(file.path(resDir, "12_allgenes_univ_surv.Rds")) == F | overwrite == T){
  print("Starting univariate analysis for overall survival")
  start <- Sys.time()
  univ_gene_formulas <- sapply(
    colnames(coxdata)[21:ncol(coxdata)],
    function(x) as.formula(
      paste(
        'Surv(time=months_overall_survival, event=overall_survival)~',
        x)))

  univ_gene_models <- lapply(univ_gene_formulas, function(x){coxph(x, data = coxdata)})
  end <- Sys.time()
  print(end-start)

  saveRDS(univ_gene_models, file=file.path(resDir, "12_allgenes_univ_surv.Rds"))
}

#### Metastasis/DRS ####

print("distant_recurrence = metastasis, including those which were metastatic at diagnosis")
print("Distant_recurrence x metastasis at diagnosis:")
print(table(colData(dds)$distant_recurrence, colData(dds)$metastasis_at_diagnosis))
#months_to_drs

#Multivariate
if(file.exists(file.path(resDir, "12_allgenes_multivar_dr.Rds")) == F | overwrite == T){
  print("Starting multivariate analysis for distant recurrence")
  start <- Sys.time()
  multivar_gene_formulas <- sapply(
    colnames(coxdata)[21:ncol(coxdata)],
    function(x) as.formula(
      paste(
        'Surv(time=months_to_drs, event=distant_recurrence)~age_at_diagnosis+year_of_diagnosis+stage+grade+PAM50+',
        x)))

  multivar_gene_models <- lapply(multivar_gene_formulas, function(x){coxph(x, data = coxdata)})
  end <- Sys.time()
  print(end-start)

  saveRDS(multivar_gene_models, file=file.path(resDir, "12_allgenes_multivar_dr.Rds"))
}

#Univariate
if(file.exists(file.path(resDir, "12_allgenes_univ_dr.Rds")) == F | overwrite == T){
  print("Starting univariate analysis for distant recurrence")
  start <- Sys.time()
  univ_gene_formulas <- sapply(
    colnames(coxdata)[21:ncol(coxdata)],
    function(x) as.formula(
      paste(
        'Surv(time=months_to_drs, event=distant_recurrence)~',
        x)))

  univ_gene_models <- lapply(univ_gene_formulas, function(x){coxph(x, data = coxdata)})
  end <- Sys.time()
  print(end-start)

  saveRDS(univ_gene_models, file=file.path(resDir, "12_allgenes_univ_dr.Rds"))
}

#### Elastic Cox ####

#A pentalized cox model to determine which groups of covariates are most predictive of survival.

#We'll want two versions: One matrix with just the gene expression, and one with gene expression plus clinical covariates.
covariates <- c("age_at_diagnosis", "year_of_diagnosis","stage","grade")
enet.gene <- coxdata[,21:ncol(coxdata)] #Gene names start at position 21
rownames(enet.gene) <- coxdata$sample_name
enet.gene <- as.matrix(enet.gene)


enet.clin <- coxdata[,c(which(colnames(coxdata) %in% covariates), 21:ncol(coxdata))]
rownames(enet.clin) <- coxdata$sample_name
enet.clin <- as.matrix(enet.clin)

#Overall survival/OS
os = data.frame(time = coxdata$months_overall_survival, event = coxdata$overall_survival)
rownames(os) = coxdata$sample_name

#os = as.matrix(os)
os = Surv(os$time, os$event)

#DRS
#When looking at drs, time point 0 (metastatic at diagnosis) is not allowed
drs = data.frame(time = coxdata$months_to_drs, event = coxdata$distant_recurrence)
rownames(drs) = coxdata$sample_name
drs = drs[!rownames(drs) %in% met_at_diag,]
#Some samples have months_to_drs = 0

drs = drs[drs$time!=0,]
#drs = as.matrix(drs)
drs = Surv(drs$time, drs$event)
head(drs)

enet.gene.drs = enet.gene[!(rownames(enet.gene) %in% met_at_diag),]
enet.clin.drs = enet.clin[!(rownames(enet.clin) %in% met_at_diag),]

nested.cox <- function(x,y,nfolds=10,nfolds_inner=10,
                       s="lambda.min",
                       alpha=0.5){

  #Based on nested.cv from the TANDEM package

  n = length(y)
  foldid = ceiling(sample(1:n)/n * nfolds)

  fitted_relative_risk = rep(NA, length(y))
  #y_hat = rep(NA, length(y))

  features = list()
  fitplots = list()
  results = list()

  for (i in 1:nfolds) {
    ind = foldid == i
    x_train = x[!ind, ]
    y_train = y[!ind]
    x_test = x[ind, ]
    n_i = length(y_train)
    foldid_i = ceiling(sample(1:n_i)/n_i * nfolds_inner)
    fit = glmnet::cv.glmnet(x_train, y_train, family = "cox", alpha = alpha,
                            foldid = foldid_i)
    fitplots[[length(fitplots) + 1]] = fit
    features[[length(features) + 1]] = coef(fit, s=s)[,1] %>% enframe() %>% filter(value!=0 & name != "(Intercept)")
    fitted_relative_risk[ind] = as.vector(glmnet::predict.cv.glmnet(fit,
                                                                    newx = x_test, s = s, alpha=alpha, type="response"))

    print(i)
  }

  results$features <- features
  results$fitplots <- fitplots

  return(results)
}


#Overall survival, genes only

if(file.exists(file.path(resDir, "12_glm_os_gene.Rds")) == F | overwrite == T){
  print("Elastic Cox for overall survival, genes only...")
  start <- Sys.time()
  glm.os <- nested.cox(x = enet.gene, y = os, s = "lambda.min", alpha = 0.5)
  end <- Sys.time()
  print(end-start)
  saveRDS(glm.os, file.path(resDir,"12_glm_os_gene.Rds"))
}

#Overall survival, including clinical covariates

if(file.exists(file.path(resDir, "12_glm_os_clin.Rds")) == F | overwrite == T){
  print("Elastic Cox for overall survival, including clinical covariates...")
  start <- Sys.time()
  glm.os.clin <- nested.cox(x = enet.clin, y = os, s = "lambda.min", alpha = 0.5)
  end <- Sys.time()
  print(end-start)
  saveRDS(glm.os.clin, file.path(resDir,"12_glm_os_clin.Rds"))
}

#DRS, genes only

if(file.exists(file.path(resDir, "12_glm_drs_gene.Rds")) == F | overwrite == T){
  print("Elastic Cox for distant recurrence (metastasis), genes only...")
  start <- Sys.time()
  glm.drs <- nested.cox(x = enet.gene.drs, y = drs, s = "lambda.min", alpha = 0.5)
  end <- Sys.time()
  print(end-start)
  saveRDS(glm.drs, file.path(resDir,"12_glm_drs_gene.Rds"))
}

#DRS, including clinical covariates
if(file.exists(file.path(resDir, "12_glm_drs_clin.Rds")) == F | overwrite == T){
  print("Elastic Cox for distant recurrence, including clinical covariates...")
  start <- Sys.time()
  glm.drs.clin <- nested.cox(x = enet.clin.drs, y = drs, s = "lambda.min", alpha = 0.5)
  end <- Sys.time()
  print(end-start)
  saveRDS(glm.drs.clin, file.path(resDir,"12_glm_drs_clin.Rds"))
}
