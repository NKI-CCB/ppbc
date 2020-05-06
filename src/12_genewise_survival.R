rm(list = ls())

library(DESeq2)
library(edgeR)
library(survival)
library(survminer)
library(RColorBrewer)
library(glmnet)
library(tidyverse)
library(here)

#Whether to overwrite existing results files

overwrite <- FALSE

set.seed(123)

#### Load data ----

#Metadata containing clinical covariates, as defined in notebook 4
survdata <- readRDS(here("data", "Rds", "04_survdata.Rds"))
#Exclude unnecessary columns, see notebook 4 for details
survdata <- survdata %>% select(-reason_death, -year_of_diagnosis)
colnames(survdata)

#Gene annotation
gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id,
                               gene_name, gene_type,
                               description = gene_description) %>% distinct()

#Gene expression
dds = readRDS(here("data/Rds/08_dds_ovr_inv_vs_rest.Rds"))

#### Prepare cox data ----

#Extract raw counts
genEx <- assay(dds)

#Subset for only those samples with enough survival data
print(paste("Samples in gene expression matrix gene expression matrix:", ncol(genEx)))
genEx <- genEx[,colnames(genEx) %in% survdata$sample_name]
print(paste("Samples with adequate metadata:", ncol(genEx)))

#TMM and log transform.
#Anecdotally, TMM does a better job than VST at dealing with outliers

dge <- DGEList(genEx, samples = survdata, group = survdata$study_group)

normTMMlog2 <- function(object){
  object = calcNormFactors(object, method="TMM")
  object = cpm(object, log=T, normalized.lib.sizes=T)
  return(object)
}

mat <- normTMMlog2(dge)

#Combine into single dataframe for Cox
stopifnot(identical(rownames(t(mat)), survdata$sample_name))
coxdata = cbind(survdata, t(mat))

head(coxdata[,1:30])

saveRDS(coxdata, here("data", "Rds", "12_coxdata.Rds"))

#### Survival functions ----

extract_genewise_cox_results <- function(cox_models, type){

  #Cox models generated as follows:
  #gene_formulas <- sapply(gene_vector,function(x) as.formula(paste('Surv(time, event)~', x)))
  #cox_models <- lapply(gene_formulas, function(x){coxph(x, data = coxdata)})

  stopifnot(type %in% c("univariate", "multivariate"))

  if(type == "univariate"){
    univ_gene_results <- lapply(cox_models,
                                function(x){
                                  f <- paste0(as.character(x$formula)[2],as.character(x$formula)[1], as.character(x$formula)[3])
                                  x <- summary(x)
                                  p.value<-signif(x$wald["pvalue"], digits=2)
                                  #wald.test<-signif(x$wald["test"], digits=2) #Unnecessary
                                  beta<-signif(x$coef[1], digits=2);#coeficient beta
                                  HR <-signif(x$coef[2], digits=2);#exp(beta)
                                  HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                                  HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                                  HR <- paste0(HR, " (",
                                               HR.confint.lower, "-", HR.confint.upper, ")")
                                  res<-c(beta, HR, p.value, f)
                                  names(res)<-c("beta", "HR (95% CI for HR)",
                                                "p.value", "call")
                                  return(res)
                                  #return(exp(cbind(coef(x),confint(x))))
                                })

    res_df <- t(as.data.frame(univ_gene_results, check.names = FALSE))
  }

  if(type == "multivariate"){
    multivar_gene_results <- lapply(cox_models,
                                    function(x){
                                      f <- paste0(as.character(x$formula)[2],as.character(x$formula)[1], as.character(x$formula)[3])
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
                                      res<-c(beta, HR, p.value, f)
                                      names(res)<-c("beta", "HR (95% CI for HR)","p.value","call")
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
    formula = paste0(formula,'age_at_diagnosis+stage+grade+surgery+radiotherapy+strata(hormonetherapy)+chemotherapy+herceptin+strata(PAM50)+')
  }

  gene_formulas <- sapply(gene_list, function(x) as.formula(paste(formula, x)))
  #return(gene_formulas)
  gene_models <- lapply(gene_formulas, function(x){coxph(x, data = data)})
  #return(gene_models)
  res <- extract_genewise_cox_results(gene_models, type=type)
  res <- res %>%
    mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
    select(gene_name, fdr, beta, `HR (95% CI for HR)`, p.value,
           ensembl_gene_id, gene_type, description, call, everything()) %>%
    arrange(fdr)
  end <- Sys.time()

  if(show_runtime){print(end-start)}
  return(res)
}

#### Cox regression ----

#Results directory
resDir = here("data", "Rds")
dir.create(resDir, showWarnings = F)
stopifnot(file.exists(resDir))

#### Overall survival ----

gene_col = which(colnames(coxdata)=="ENSG00000000003")
print(paste("Gene columns from", gene_col, "to", ncol(coxdata), "of coxdata"))
print(head(coxdata[,1:30]))

#Multivariate

resPath = file.path(resDir, "12_multi_genewise_os.Rds")

#Test before running
#testdata = coxdata[1:80, 1:30]
#genewise_cox(
#  gene_list = colnames(testdata)[gene_col:ncol(testdata)],
#  data = testdata,
#  show_runtime = T,
#  time="months_overall_survival",
#  event="overall_survival",
#  type="multivariate")


if(file.exists(resPath) == F | overwrite == T){
  print("Multivariate survival")
  res <- 
  genewise_cox(
    gene_list = colnames(coxdata)[gene_col:ncol(coxdata)], #colnames(testdata)[gene_col:ncol(testdata)]
    data = coxdata, #testdata
    show_runtime = T,
    time="months_overall_survival",
    event="overall_survival",
    type="multivariate")
  saveRDS(res, resPath)
}

#Univariate
resPath = file.path(resDir, "12_uni_genewise_os.Rds")

#Test before running
#genewise_cox(
#  gene_list = colnames(testdata)[gene_col:ncol(testdata)],
#  data = testdata,
#  show_runtime = T,
#  time="months_overall_survival",
#  event="overall_survival",
#  type="univariate")

if(file.exists(resPath) == F| overwrite == T){
  print("Univariate survival")
  res <- genewise_cox(
    gene_list = colnames(coxdata)[gene_col:ncol(coxdata)],
    data = coxdata,
    show_runtime = T,
    time="months_overall_survival",
    event="overall_survival",
    type="univariate")
  saveRDS(res, resPath)
}

#### Metastasis/DRS ####

#print("distant_recurrence = metastasis, including those which were metastatic at diagnosis")
#print("Distant_recurrence x metastasis at diagnosis:")
#print(table(colData(dds)$distant_recurrence, colData(dds)$metastasis_at_diagnosis))
#months_to_drs

#Multivariate
resPath = file.path(resDir, "12_multi_genewise_drs.Rds")

#Testing
#genewise_cox(
#  gene_list = colnames(testdata)[gene_col:ncol(testdata)],
#  data = testdata,
#  show_runtime = T,
#  time="months_to_drs",
#  event="distant_recurrence",
#  type="multivariate")

if(file.exists(resPath) == F| overwrite == T){
  print("Multivariate metastasis")
  res <- genewise_cox(
    gene_list = colnames(coxdata)[gene_col:ncol(coxdata)],
    data = coxdata,
    show_runtime = T,
    time="months_to_drs",
    event="distant_recurrence",
    type="multivariate")
  saveRDS(res, resPath)
}

#Univariate
resPath = file.path(resDir, "12_uni_genewise_drs.Rds")

#Testing
#genewise_cox(
#  gene_list = colnames(testdata)[gene_col:ncol(testdata)],
#  data = testdata,
#  show_runtime = T,
#  time="months_to_drs",
#  event="distant_recurrence",
#  type="univariate")

if(file.exists(resPath) == F| overwrite == T){
  print("Univariate metastasis")
  res <- genewise_cox(
    gene_list = colnames(coxdata)[gene_col:ncol(coxdata)],
    data = coxdata,
    show_runtime = T,
    time="months_to_drs",
    event="distant_recurrence",
    type="univariate")
  saveRDS(res, resPath)
}

#### Involution samples only ----

#Subset for only those samples with enough survival data
print("All samples")
print(ncol(genEx))

print("Involution samples:")
invEx <- genEx[,colnames(genEx)[str_detect(colnames(genEx), "inv")]]
ncol(invEx)

invmeta <- survdata %>% filter(sample_name %in% colnames(invEx))

#TMM and log transform

invdge <- DGEList(invEx, samples = invmeta, group = invmeta$study_group)
invmat <- normTMMlog2(invdge)

#Combine into single dataframe for Cox
stopifnot(identical(rownames(t(invmat)), invmeta$sample_name))

invdata = cbind(invmeta, t(invmat))
inv_gene_col = which(colnames(invdata)=="ENSG00000000003")

print(paste("Gene columns from", inv_gene_col, "to", ncol(invdata), "of involution-only coxdata"))
print(head(invdata[,1:30]))
saveRDS(invdata, file.path(resDir, "12_invdata.Rds"))

#### Overall survival ----

#Testing
#testdata = invdata[,1:30]
#genewise_cox(
#  gene_list = colnames(testdata)[inv_gene_col:ncol(testdata)],
#  data = testdata,
#  show_runtime = T,
#  time="months_overall_survival",
#  event="overall_survival",
#  type="multivariate")

resPath = file.path(resDir, "12_inv_multi_genewise_os.Rds")

if(file.exists(resPath) == F| overwrite == T){
  print("Multivariate survival, involution samples")
  res <- genewise_cox(
    gene_list = colnames(invdata)[inv_gene_col:ncol(invdata)],
    data = invdata,
    show_runtime = T,
    time="months_overall_survival",
    event="overall_survival",
    type="multivariate")
  saveRDS(res, resPath)
}

#Univariate
resPath = file.path(resDir, "12_inv_uni_genewise_os.Rds")

#Testing
#testdata = invdata[,1:30]
#genewise_cox(
#  gene_list = colnames(testdata)[inv_gene_col:ncol(testdata)],
#  data = testdata,
#  show_runtime = T,
#  time="months_overall_survival",
#  event="overall_survival",
#  type="univariate")

if(file.exists(resPath) == F| overwrite == T){
  print("Univariate survival, involution samples")
  res <- genewise_cox(
    gene_list = colnames(invdata)[inv_gene_col:ncol(invdata)],
    data = invdata,
    show_runtime = T,
    time="months_overall_survival",
    event="overall_survival",
    type="univariate")
  saveRDS(res, resPath)
}

#### DRS ----

#Multivariate
resPath = file.path(resDir, "12_inv_multi_genewise_drs.Rds")

#genewise_cox(
#  gene_list = colnames(testdata)[inv_gene_col:ncol(testdata)],
#  data = testdata,
#  show_runtime = T,
#  time="months_to_drs",
#  event="distant_recurrence",
#  type="multivariate")

if(file.exists(resPath) == F| overwrite == T){
  print("Multivariate metastasis, involution samples")
  res <- genewise_cox(
    gene_list = colnames(invdata)[inv_gene_col:ncol(invdata)],
    data = invdata,
    show_runtime = T,
    time="months_to_drs",
    event="distant_recurrence",
    type="multivariate")
  saveRDS(res, resPath)
}

#Univariate
resPath = file.path(resDir, "12_inv_uni_genewise_drs.Rds")

#genewise_cox(
#  gene_list = colnames(testdata)[inv_gene_col:ncol(testdata)],
#  data = testdata,
#  show_runtime = T,
#  time="months_to_drs",
#  event="distant_recurrence",
#  type="univariate")

if(file.exists(resPath) == F| overwrite == T){
  print("Univariate metastasis, involution samples")
  res <- genewise_cox(
    gene_list = colnames(invdata)[inv_gene_col:ncol(invdata)],
    data = invdata,
    show_runtime = T,
    time="months_to_drs",
    event="distant_recurrence",
    type="univariate")
  saveRDS(res, resPath)
}

