library(survival)
library(survminer)
library(RColorBrewer)
library(glmnet)
library(tidyverse)
library(here)

#Whether to overwrite existing results files

overwrite <- T

set.seed(123)

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

# Allow the file to be sourced for the functions above,
# without running the code below
if(sys.nframe() == 0){
  
  #Load data
  coxdata <- readRDS(here("data/rnaseq/processed/12_coxdata.Rds"))
  gx_annot <- read_tsv(here("data/rnaseq/metadata/01_gene_annot.tsv")) %>% 
    select(ensembl_gene_id = gene_id,
           gene_name, gene_type,
           description = gene_description) %>% distinct()
  
  
  #Results directory
  resDir = here("data/rnaseq/processed")
  dir.create(resDir, showWarnings = F)
  stopifnot(file.exists(resDir))
  
  #### Overall survival ----
  
  gene_col = which(colnames(coxdata)=="ENSG00000000003")
  print(paste("Gene columns from", gene_col, "to", ncol(coxdata), "of coxdata: all samples"))
  print(head(coxdata[,(gene_col-3):(gene_col+1)]))
  
  #Multivariate
  
  resPath = file.path(resDir, "12_multi_genewise_os.Rds")
  
  # Test before running
  # testdata = coxdata[1:80, 1:30]
  # genewise_cox(
  #  gene_list = colnames(testdata)[gene_col:ncol(testdata)],
  #  data = testdata,
  #  show_runtime = T,
  #  time="time_OS_months",
  #  event="death",
  #  type="multivariate")
  
  
  if(file.exists(resPath) == F | overwrite == T){
    print("Multivariate survival: all samples")
    res <- 
      genewise_cox(
        gene_list = colnames(coxdata)[gene_col:ncol(coxdata)], #colnames(testdata)[gene_col:ncol(testdata)]
        data = coxdata, #testdata
        show_runtime = T,
        time="time_OS_months",
        event="death",
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
  #  time="time_OS_months",
  #  event="death",
  #  type="univariate")
  
  if(file.exists(resPath) == F | overwrite == T){
    print("Univariate survival: all samples")
    res <- genewise_cox(
      gene_list = colnames(coxdata)[gene_col:ncol(coxdata)],
      data = coxdata,
      show_runtime = T,
      time="time_OS_months",
      event="death",
      type="univariate")
    saveRDS(res, resPath)
  }
  
  #### Metastasis/DRS ####
  
  #print("distant_recurrence = metastasis, including those which were metastatic at diagnosis")
  #print("Distant_recurrence x metastasis at diagnosis:")
  #print(table(colData(dds)$distant_recurrence, colData(dds)$metastasis_at_diagnosis))
  #time_DRS_months
  
  #Multivariate
  resPath = file.path(resDir, "12_multi_genewise_drs.Rds")
  
  #Testing
  #genewise_cox(
  #  gene_list = colnames(testdata)[gene_col:ncol(testdata)],
  #  data = testdata,
  #  show_runtime = T,
  #  time="time_DRS_months",
  #  event="distant_recurrence",
  #  type="multivariate")
  
  if(file.exists(resPath) == F | overwrite == T){
    print("Multivariate metastasis: all samples")
    res <- genewise_cox(
      gene_list = colnames(coxdata)[gene_col:ncol(coxdata)],
      data = coxdata,
      show_runtime = T,
      time="time_DRS_months",
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
  #  time="time_DRS_months",
  #  event="distant_recurrence",
  #  type="univariate")
  
  if(file.exists(resPath) == F | overwrite == T){
    print("Univariate metastasis: all samples")
    res <- genewise_cox(
      gene_list = colnames(coxdata)[gene_col:ncol(coxdata)],
      data = coxdata,
      show_runtime = T,
      time="time_DRS_months",
      event="distant_recurrence",
      type="univariate")
    saveRDS(res, resPath)
  }
  
}
