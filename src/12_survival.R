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

#Created in notebook 12_survival.Rmd
coxdata <- readRDS(here("data", "Rds", "12_coxdata.Rds"))
#Gene annotation
gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>% distinct()
dds = readRDS(here("data/Rds/08_dds_inv_vs_rest_standard.Rds"))


extract_genewise_cox_results <- function(cox_models, type){

  #Cox models generated as follows:
  #gene_formulas <- sapply(gene_vector,function(x) as.formula(paste('Surv(time, event)~', x)))
  #cox_models <- lapply(gene_formulas, function(x){coxph(x, data = coxdata)})

  stopifnot(type %in% c("univariate", "multivariate"))

  if(type == "univariate"){
    univ_gene_results <- lapply(cox_models,
                                function(x){
                                  x <- summary(x)
                                  p.value<-signif(x$wald["pvalue"], digits=2)
                                  #wald.test<-signif(x$wald["test"], digits=2) #Unnecessary
                                  beta<-signif(x$coef[1], digits=2);#coeficient beta
                                  HR <-signif(x$coef[2], digits=2);#exp(beta)
                                  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                                  HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                                  HR <- paste0(HR, " (",
                                               HR.confint.lower, "-", HR.confint.upper, ")")
                                  res<-c(beta, HR,
                                         #wald.test,
                                         p.value)
                                  names(res)<-c("beta", "HR (95% CI for HR)", #"wald.test",
                                                "p.value")
                                  return(res)
                                  #return(exp(cbind(coef(x),confint(x))))
                                })

    res_df <- t(as.data.frame(univ_gene_results, check.names = FALSE))
  }

  if(type == "multivariate"){
    multivar_gene_results <- lapply(cox_models,
                                    function(x){
                                      x <- summary(x)
                                      #Drop the other formula elements
                                      coefs <- as.data.frame(x$coefficients)
                                      coefs <- rownames_to_column(coefs, "ensembl_gene_id")
                                      gene = tail(coefs$ensembl_gene_id,1)
                                      coefs <- coefs[coefs$ensembl_gene_id==gene, , drop=F]
                                      colnames(coefs) <- c("ensembl_gene_id", "beta", "HR", "se(beta)", "z", "p.value")

                                      #Extract p val and stats for gene of interest
                                      p.value<-signif(coefs$p.value, digits = 2)
                                      beta<-signif(coefs$beta, digits=2);#coeficient beta
                                      HR <-signif(coefs$HR, digits=2);#exp(beta)

                                      #Do the same for conf intervals
                                      conf.int = rownames_to_column(as.data.frame(x$conf.int), "ensembl_gene_id")
                                      conf.int = conf.int[conf.int$ensembl_gene_id==gene, , drop=F]
                                      HR.confint.lower <- signif(conf.int[,"lower .95"], 2)
                                      HR.confint.upper <- signif(conf.int[,"upper .95"],2)
                                      HR <- paste0(HR, " (",
                                                   HR.confint.lower, "-", HR.confint.upper, ")")
                                      res<-c(beta, HR, p.value)
                                      names(res)<-c("beta", "HR (95% CI for HR)","p.value")
                                      return(res)
                                    })

    res_df <- t(as.data.frame(multivar_gene_results, check.names = FALSE))
  }

  df = as.data.frame(res_df) %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(., gx_annot,by="ensembl_gene_id") %>%
    mutate(p.value = as.numeric(as.character(p.value)),
           beta = as.numeric(as.character(beta))) %>%
    select(gene_name, p.value, everything()) %>%
    arrange(p.value)

  return(df)
}

genewise_cox <- function(gene_list, time, event, data = coxdata, type, show_runtime=T){
  stopifnot(type %in% c("univariate", "multivariate"))
  stopifnot(time %in% colnames(data))
  stopifnot(event %in% colnames(data))

  start <- Sys.time()

  formula = paste0('Surv(time=',time,', event=',event,')~')

  if(type == "multivariate"){
    formula = paste0(formula,'age_at_diagnosis+year_of_diagnosis+stage+grade+surgery+radiotherapy+hormonetherapy+chemotherapy+herceptin+PAM50+')
  }

  gene_formulas <- sapply(gene_list, function(x) as.formula(paste(formula, x)))
  gene_models <- lapply(gene_formulas, function(x){coxph(x, data = data)})
  res <- extract_genewise_cox_results(gene_models, type=type)
  end <- Sys.time()

  if(show_runtime){print(end-start)}
  return(res)
}

#### Cox regression ####

#Results directory
resDir = here("data", "Rds")
dir.create(resDir, showWarnings = F)
stopifnot(file.exists(resDir))

#### Overall survival ####

#Multivariate

#Old version
#if(file.exists(file.path(resDir, "12_allgenes_multivar_surv.Rds")) == F | overwrite == T){
#  print("Starting multivariate analysis for overall survival")
#  start <- Sys.time()
#  multivar_gene_formulas <- sapply(
#    colnames(coxdata)[21:ncol(coxdata)],
#    function(x) as.formula(
#      paste(
#        'Surv(time=months_overall_survival, event=overall_survival)~age_at_diagnosis+year_of_diagnosis+stage+grade+PAM50+',
#        x)))

#  multivar_gene_models <- lapply(multivar_gene_formulas, function(x){coxph(x, data = coxdata)})
#  end <- Sys.time()
#  print(end-start)

#  saveRDS(multivar_gene_models, file=file.path(resDir, "12_allgenes_multivar_surv.Rds"))
#}

resPath = file.path(resDir, "12_allgenes_multivar_surv.Rds")

if(file.exists(resPath) == F| overwrite == T){
  print("Multivariate survival")
  res <- genewise_cox(
    gene_list = colnames(coxdata)[21:ncol(coxdata)],
    data = coxdata,
    show_runtime = T,
    time="months_overall_survival",
    event="overall_survival",
    type="multivariate")
  saveRDS(res, resPath)
}

#Univariate
resPath = file.path(resDir, "12_allgenes_univ_surv.Rds")

if(file.exists(resPath) == F| overwrite == T){
  print("Univariate survival")
  res <- genewise_cox(
    gene_list = colnames(coxdata)[21:ncol(coxdata)],
    data = coxdata,
    show_runtime = T,
    time="months_overall_survival",
    event="overall_survival",
    type="univariate")
  saveRDS(res, resPath)
}

#### Metastasis/DRS ####

print("distant_recurrence = metastasis, including those which were metastatic at diagnosis")
print("Distant_recurrence x metastasis at diagnosis:")
print(table(colData(dds)$distant_recurrence, colData(dds)$metastasis_at_diagnosis))
#months_to_drs

#Multivariate
resPath = file.path(resDir, "12_allgenes_multivar_dr.Rds")

if(file.exists(resPath) == F| overwrite == T){
  print("Multivariate metastasis")
  res <- genewise_cox(
    gene_list = colnames(coxdata)[21:ncol(coxdata)],
    data = coxdata,
    show_runtime = T,
    time="months_to_drs",
    event="distant_recurrence",
    type="multivariate")
  saveRDS(res, resPath)
}

#Univariate
resPath = file.path(resDir, "12_allgenes_univ_dr.Rds")

if(file.exists(resPath) == F| overwrite == T){
  print("Univariate metastasis")
  res <- genewise_cox(
    gene_list = colnames(coxdata)[21:ncol(coxdata)],
    data = coxdata,
    show_runtime = T,
    time="months_to_drs",
    event="distant_recurrence",
    type="univariate")
  saveRDS(res, resPath)
}


