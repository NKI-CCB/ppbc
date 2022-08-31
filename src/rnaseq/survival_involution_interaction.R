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
overwrite <- T

#### Cox functions ----

extract_interaction_cox_results <- function(cox_models, anno_df = gx_annot){
  
  gene_results <- lapply(cox_models,
                         function(x){
                           
                           #Save formula
                           f <- paste0(as.character(x$formula)[2],as.character(x$formula)[1], as.character(x$formula)[3])
                           
                           x <- summary(x)
                           
                           #Overall model significance
                           logrank_p <- x$logtest[[3]] #p value logrank test
                           #wald_p <- x$waldtest[[3]] #p value wald test
                           #lrt_p <- x$sctest[[3]] #p value likelihood ratio test
                           
                           #Coefficients
                           coefs <- as.data.frame(x$coefficients)
  
                           #Save gene from coefficient
                           coefs <- rownames_to_column(coefs, "feature")
                           gene <- coefs$feature[str_detect(coefs$feature, "^ENSG")]
                           
                           #Extract results for interaction coefficient
                           
                           int <- coefs[coefs$feature==coefs$feature[str_detect(coefs$feature, ":")], , drop=F]
                           colnames(int) <- c("feature", "beta", "HR", "se(beta)", "z", "p.value")
                           int <- int %>% mutate(p.value = signif(p.value, 2),
                                                 beta = signif(beta, 2),
                                                 HR = signif(HR, 2),
                                                 ensembl_gene_id = gene,
                                                 interaction = feature,
                                                 logrank_p = signif(logrank_p, 2),
                                                 #wald_p = signif(logrank_p, 2),
                                                 #lrt_p = signif(logrank_p, 2),
                                                 call = f
                           )
                           
                           #Do the same for conf intervals
                           conf = rownames_to_column(as.data.frame(x$conf.int), "feature")
                           conf = conf[conf$feature == coefs$feature[str_detect(coefs$feature, ":")], , drop=F]
                           HR.confint.lower <- signif(conf[,"lower .95"], 2)
                           HR.confint.upper <- signif(conf[,"upper .95"],2)
                           int$HR <- paste0(int$HR, " (", 
                                            HR.confint.lower, "-", HR.confint.upper, ")")
                           res <- int %>% select(ensembl_gene_id, beta, HR, p.value, everything())
                           res <- res %>% dplyr::rename(`HR (95% CI for HR)` = HR)
                           res <- res %>% select(-feature, -`se(beta)`, -z) #redundant
                           res <- as.vector(res[1,,drop=T])
                           
                           return(res)
                         })
  #res_df <- t(as.data.frame(gene_results, check.names = FALSE))
  #return(gene_results)
  res_df <- bind_rows(gene_results)
  #return(res_df)
  
  df = res_df %>%
    left_join(., anno_df, by="ensembl_gene_id") %>%
    mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
    select(gene_name, fdr, beta,`HR (95% CI for HR)`, p.value, gene_type, description, everything())
  
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
  
  res <- extract_interaction_cox_results(gene_models)
  
  #We added this already
  #res$formula <- gsub("[[:space:]]", "", unlist(format(gene_formulas)))
  
  res <- res %>%arrange(fdr, p.value)
  
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
  
  gene_col = which(colnames(coxdata)=="ENSG00000000003")
  print(paste("Gene columns from", gene_col, "to", ncol(coxdata), "of coxdata"))
  print(head(coxdata[,(gene_col-3):(gene_col+1)]))
  
  #### Test dataset ----
  
  # #First test multivariate models
  # test = genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:(gene_col + 6)],
  #                                  survival_type = "os", formula_type = "multivariate",
  #                                  return_models=T)
  # test[[1]]
  # test2 = extract_interaction_cox_results(test)
  # test2$call[1]
  # 
  # #When return_models is false, extract_interaction_cox_results() is called as part of genewise_cox_interactions()
  # test3 = genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:(gene_col + 6)],
  #                                   survival_type = "os", formula_type = "multivariate",
  #                                   return_models=F)
  # test3
  # 
  # #Then univariate
  # test4 = genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:(gene_col + 6)],
  #                                   survival_type = "os", formula_type = "univariate",
  #                                   return_models=T)
  # test4[[1]]
  # 
  # test5 = extract_interaction_cox_results(test4)
  # test5
  # 
  # #When return_models is false, extract_interaction_cox_results() is called as part of genewise_cox_interactions()
  # test6=genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:(gene_col + 6)],
  #                                 survival_type = "os", formula_type = "univariate",
  #                                 return_models=F)
  # test6
  # 
  # rm(list=c("test","test2","test3","test4","test5","test6"))
  # 
  
  #### Cox regression ----
  
  #Results directory
  resDir = here("data", "rnaseq", "processed")
  dir.create(resDir, showWarnings = F)
  stopifnot(file.exists(resDir))
  
  #### OS ####
  
  #Univariate
  resPath = file.path(resDir, "13_uni_interaction_os.Rds")
  
  if(file.exists(resPath) == F| overwrite == T){
    print("Univariate interaction involution*gene for overall survival")
    res <-  genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:ncol(coxdata)], 
                                      survival_type = "os",
                                      formula_type = "univariate")
    saveRDS(res, resPath)
  }
  
  #Multivariate
  resPath = file.path(resDir, "13_multi_interaction_os.Rds")
  
  if(file.exists(resPath) == F| overwrite == T){
    print("Multivariate interaction involution*gene for overall survival")
    res <-  genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:ncol(coxdata)], 
                                      survival_type = "os",
                                      formula_type = "multivariate")
    saveRDS(res, resPath)
  }
  
  #### DRS ####
  
  #Univariate
  resPath = file.path(resDir, "13_uni_interaction_drs.Rds")
  
  if(file.exists(resPath) == F| overwrite == T){
    print("Univariate interaction involution*gene for DRS")
    res <-  genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:ncol(coxdata)], 
                                      survival_type = "drs",
                                      formula_type = "univariate")
    saveRDS(res, resPath)
  }
  
  #Multivariate
  resPath = file.path(resDir, "13_multi_interaction_drs.Rds")
  
  if(file.exists(resPath) == F| overwrite == T){
    print("Multivariate interaction involution:gene for DRS")
    res <-  genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:ncol(coxdata)], 
                                      survival_type = "drs",
                                      formula_type = "multivariate")
    saveRDS(res, resPath)
  }
}