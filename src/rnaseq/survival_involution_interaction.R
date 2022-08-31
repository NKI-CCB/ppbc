library(DESeq2)
library(survival)
library(survminer)
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(glmnet)
library(openxlsx)
library(here)
library(ggpubr)
library(tidyverse)

# Whether to overwrite pre-existing results
overwrite <- F

#### Cox functions ----

tidy_cox <- function(cox_models, anno_df = gx_annot){
  
  #Extract results
  res_df <- map_df(cox_models, broom::tidy, conf.int = T, .id = "ensembl_gene_id") %>%
    filter(str_detect(term, "involution:ENSG")) %>%
    mutate(beta = signif(estimate, 2),
           HR = paste0(signif(exp(estimate), 2),
                       " (", signif(conf.low, 2), "-",
                       signif(conf.high,2), ")")) %>%
    select(ensembl_gene_id, p.value, beta, HR, interaction = term)
  
  calls <- sapply(cox_models, function(x){
    f <- paste(as.character(x$formula)[2],as.character(x$formula)[1], as.character(x$formula)[3])
  }) %>% enframe("ensembl_gene_id", "call")
  
  res_df <- left_join(res_df, calls, by = "ensembl_gene_id")
  
  df = res_df %>%
    left_join(., anno_df, by="ensembl_gene_id") %>%
    mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
    select(gene_name, fdr, beta,`HR (95% CI for HR)`=HR, p.value,
           gene_type, description, everything())
  
  return(df)
}

genewise_cox_interactions <- function(gene_list, data = coxdata, survival_type,
                                      formula_type, show_runtime=T, return_models=F){
  
  stopifnot(survival_type %in% c("os", "drs"))
  stopifnot(formula_type %in% c("univariate", "multivariate"))
  
  start <- Sys.time()
  
  if(survival_type == "os"){
    time = "time_OS_months"
    event = "death"
  } else {
    time = "time_DRS_months"
    event = "distant_recurrence"
  }
  
  formula = paste0('Surv(time=',time,', event=',event,')~')
  
  if (formula_type == "multivariate"){
    cov = 'age_at_diagnosis+stage+grade+surgery+radiotherapy+strata(hormonetherapy)+chemotherapy+herceptin+strata(PAM50)'
    formula = paste0(formula, cov)
    gene_formulas <-  sapply(gene_list,
                             function(x)
                               as.formula(paste(formula, paste0("involution*", x), sep="+")))
  } else {
    gene_formulas <-  sapply(gene_list,
                             function(x)
                               as.formula(paste(formula, paste0("involution*", x))))
  }
  
  if(show_runtime){print(formula)}
  
  #return(unlist(format(gene_formulas)))
  gene_models <- lapply(gene_formulas, function(x){coxph(x, data = data)})
  
  if(return_models==T){
    return(gene_models)  
  }
  
  res <- tidy_cox(gene_models)
  
  #We added this already
  #res$formula <- gsub("[[:space:]]", "", unlist(format(gene_formulas)))
  
  res <- res %>% arrange(fdr, p.value)
  
  end <- Sys.time()
  
  if(show_runtime){print(end-start)}
  
  return(res)
}

# Allow the file to be sourced for the functions above,
# without running the code below
if(sys.nframe() == 0){

  #### Load data ----
  
  #Gene annotation
  gx_annot <- read_tsv(here("data/rnaseq/metadata/01_gene_annot.tsv"))
  gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id,
                                 gene_name, gene_type,
                                 description = gene_description) %>% distinct()
  
  #Clinical covariates and TMM-normalized/log transformed gene expression
  #In samples x features format
  coxdata <- readRDS(here("data/rnaseq/processed/12_coxdata.Rds"))
  colnames(coxdata)[1:30]
  
  ## Involuting vs rest encouding
  levels(coxdata$study_group)
  coxdata = coxdata %>%
    mutate(involution = as.numeric(study_group == "ppbcpw", .after = study_group))
  
  table(coxdata$study_group, coxdata$involution)
  
  #Detect start of count data
  gene_col = which(colnames(coxdata)=="ENSG00000000003")
  print(paste("Gene columns from", gene_col, "to", ncol(coxdata), "of coxdata"))
  print(head(coxdata[,(gene_col-3):(gene_col+1)]))
  
  #### Test dataset ----
  
  #First test multivariate models
  # test = genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:(gene_col + 6)],
  #                                  survival_type = "os", formula_type = "multivariate",
  #                                  return_models=T)
  # test
  # tidy_cox(test) %>% arrange(fdr, p.value)
  # 
  # 
  # #When return_models is false, tidy_cox() is called as part of genewise_cox_interactions()
  # test2 = genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:(gene_col + 6)],
  #                                   survival_type = "os", formula_type = "multivariate",
  #                                   return_models=F)
  # 
  # test2
  # 
  # #Then univariate
  # test3 = genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:(gene_col + 6)],
  #                                   survival_type = "os", formula_type = "univariate",
  #                                   return_models=T)
  # test3
  # 
  # tidy_cox(test3) %>% arrange(fdr, p.value)
  # 
  # test4=genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:(gene_col + 6)],
  #                                 survival_type = "os", formula_type = "univariate",
  #                                 return_models=F)
  # test4
  # 
  # rm(list=c("test","test2","test3","test4"))

  #### Cox regression ----
  
  #Results directory
  resDir = here("data", "rnaseq", "processed")
  dir.create(resDir, showWarnings = F)
  stopifnot(file.exists(resDir))
  
  #### OS ####
  
  #Univariate
  resPath = file.path(resDir, "12_uni_interaction_os.Rds")
  
  if(file.exists(resPath) == F | overwrite == T){
    print("Univariate interaction involution*gene for overall survival")
    res <-  genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:ncol(coxdata)],
                                      survival_type = "os",
                                      formula_type = "univariate")
    saveRDS(res, resPath)
  }
  
  #Multivariate
  resPath = file.path(resDir, "12_multi_interaction_os.Rds")
  
  if(file.exists(resPath) == F| overwrite == T){
    print("Multivariate interaction involution*gene for overall survival")
    res <-  genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:ncol(coxdata)], 
                                      survival_type = "os",
                                      formula_type = "multivariate")
    saveRDS(res, resPath)
  }
  
  #### DRS ####
  
  #Univariate
  resPath = file.path(resDir, "12_uni_interaction_drs.Rds")
  
  if(file.exists(resPath) == F| overwrite == T){
    print("Univariate interaction involution*gene for DRS")
    res <-  genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:ncol(coxdata)], 
                                      survival_type = "drs",
                                      formula_type = "univariate")
    saveRDS(res, resPath)
  }
  
  #Multivariate
  resPath = file.path(resDir, "12_multi_interaction_drs.Rds")
  
  if(file.exists(resPath) == F| overwrite == T){
    print("Multivariate interaction involution:gene for DRS")
    res <-  genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:ncol(coxdata)], 
                                      survival_type = "drs",
                                      formula_type = "multivariate")
    saveRDS(res, resPath)
  }
}

