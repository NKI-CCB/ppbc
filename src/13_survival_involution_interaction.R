
rm(list = ls())

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

#### Load data ----

#Gene annotation
gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id,
                               gene_name, gene_type,
                               description = gene_description) %>% distinct()

#Clinical covariates and TMM-normalized/log transformed gene expression
#In samples x features format
coxdata <- readRDS(here("data", "Rds", "12_coxdata.Rds"))
colnames(coxdata)[1:30]

## Involuting vs rest encouding
levels(coxdata$study_group)
coxdata = coxdata %>%
  mutate(involution = as.numeric(study_group == "ppbc_inv")) %>%
  select(sample_name:study_group, involution, everything())

table(coxdata$study_group, coxdata$involution)

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
                           
                           int <- coefs[coefs$feature==paste0("involution:", gene), , drop=F]
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
                           conf = conf[conf$feature == paste0("involution:", gene), , drop=F]
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

genewise_cox_interactions <- function(gene_list, data = coxdata, survival_type, show_runtime=T){
  
  stopifnot(survival_type %in% c("os", "drs"))
  
  start <- Sys.time()
  
  if(survival_type == "os"){
    time = "months_overall_survival"
    event = "overall_survival"
  } else {
    time = "months_to_drs"
    event = "distant_recurrence"
  }
  
  cov = 'age_at_diagnosis+stage+grade+surgery+radiotherapy+strata(hormonetherapy)+chemotherapy+herceptin+strata(PAM50)+involution'
  
  formula = paste0('Surv(time=',time,', event=',event,')~',cov)
  
  if(show_runtime){print(formula)}
  
  gene_formulas <-  sapply(gene_list,
                           function(x)
                             as.formula(paste(formula, x, paste0("involution*", x), sep="+")))
  #return(unlist(format(gene_formulas)))
  gene_models <- lapply(gene_formulas, function(x){coxph(x, data = data)})
  #return(gene_models)
  
  res <- extract_interaction_cox_results(gene_models)
  
  #We added this already
  #res$formula <- gsub("[[:space:]]", "", unlist(format(gene_formulas)))
  
  res <- res %>%arrange(fdr, p.value)
  
  end <- Sys.time()
  
  if(show_runtime){print(end-start)}
  
  return(res)
}


#### Test dataset --- 

gene_col = which(colnames(coxdata)=="ENSG00000000003")
print(paste("Gene columns from", gene_col, "to", ncol(coxdata), "of coxdata"))

#head(coxdata[,(gene_col-2):(gene_col+2)])

test = genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:(gene_col + 6)],
                                 survival_type = "os")

test

#### Cox regression ----

#Results directory
resDir = here("data", "Rds")
dir.create(resDir, showWarnings = F)
stopifnot(file.exists(resDir))

#### OS ####
resPath = file.path(resDir, "13_multi_interaction_os.Rds")

if(file.exists(resPath) == F| overwrite == T){
  print("Interaction involution:gene for overall survival")
  res <-  genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:ncol(coxdata)], 
                                    survival_type = "os")
  saveRDS(res, resPath)
}

#### DRS ####

resPath = file.path(resDir, "13_multi_interaction_drs.Rds")

if(file.exists(resPath) == F| overwrite == T){
  print("Interaction involution:gene for DRS")
  res <-  genewise_cox_interactions(gene_list = colnames(coxdata)[gene_col:ncol(coxdata)], 
                                    survival_type = "drs")
  saveRDS(res, resPath)
}
