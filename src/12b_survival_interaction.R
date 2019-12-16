
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

# Pre-process and load data

source(here("src", "deseq_report_functions.R"))

## Count and sample data
dds = readRDS(here("data/Rds/08_dds_inv_vs_rest_standard.Rds"))
gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>% distinct()



vsd = readRDS(here("data","Rds","08_vsd_standard_vsrest.Rds"))

vsd = rownames_to_column(as.data.frame(assay(vsd)), "ensembl_gene_id")
vsd = right_join(select(gx_annot, gene_name, ensembl_gene_id), vsd, by = "ensembl_gene_id") %>%
  select(-ensembl_gene_id)
vsd = summarize_expression_duplicate_ids(vsd, id_column = "gene_name", verbose = T)

## Gene sets

#For Fisher's exact tests
gene_sets = list(
  go_bp = flexgsea::read_gmt(here("data","external","gmt","c5.bp.v7.0.symbols.gmt")),
  hallmark = flexgsea::read_gmt(here("data","external","gmt","h.all.v7.0.symbols.gmt")),
  c2_canon = flexgsea::read_gmt(here("data","external","gmt","c2.cp.v7.0.symbols.gmt")),
  c2_cgp = flexgsea::read_gmt(here("data","external","gmt","c2.cgp.v7.0.symbols.gmt"))
)


## Heatmap colors

#Created in notebook 6
all_colors = readRDS(here("data","Rds", "06_heatmap_colors.Rds"))
study_colors = all_colors$study_colors
pam_colors = all_colors$pam_colors
gene_colors = all_colors$gene_colors

## Survival metadata

#As described in notebook 12, we have 183 samples with adequate data.

metadata = readxl::read_excel(here("data", "metadata", "12_survival_metadata.xlsx"))
nrow(metadata)

## Transposed matrix for Cox regressions
#Also created in notebook 12.

coxdata = readRDS(here("data", "Rds", "12_coxdata.Rds"))
colnames(coxdata)[1:30]


## Separate columns per PPBC level
levels(coxdata$study_group)

coxdata = coxdata %>%
  mutate(non_prbc = as.numeric(study_group == "non_prbc"),
         prbc = as.numeric(study_group == "prbc"),
         ppbc_lac = as.numeric(study_group == "ppbc_lac"),
         ppbc_inv = as.numeric(study_group == "ppbc_inv")) %>%
  select(sample_name:study_group, non_prbc:ppbc_inv, everything())

#head(coxdata)

#colnames(coxdata)[1:30]


# Functions for interaction term Cox regressions

## Extract results from models
extract_interaction_cox_results <- function(cox_models, type = "full", interaction_coef = "ppbc_inv",
                                            anno_df = gx_annot, time, event){

  stopifnot(type == "full") #"interaction_only" not yet implemented


  if(type == "full"){
    gene_results <- lapply(cox_models,
                           function(x){
                             x <- summary(x)

                             #Overall model significance
                             logrank_p <- x$logtest[[3]] #p value logrank test
                             wald_p <- x$waldtest[[3]] #p value wald test
                             lrt_p <- x$sctest[[3]] #p value likelihood ratio test

                             #Save formula
                             coefs <- as.data.frame(x$coefficients)
                             form = paste0(paste0('Surv(time=',time,', event=',event,')~'), paste0(rownames(coefs), collapse = "+"))

                             #Save gene from formula
                             coefs <- rownames_to_column(coefs, "feature")
                             gene = coefs$feature[str_detect(coefs$feature, "^ENSG")]

                             #Extract results for interaction coefficient

                             int <- coefs[coefs$feature==paste0(interaction_coef, ":", gene), , drop=F]
                             colnames(int) <- c("feature", "beta", "HR", "se(beta)", "z", "p.value")
                             int <- int %>% mutate(p.value = signif(p.value, 2),
                                                   beta = signif(beta, 2),
                                                   HR = signif(HR, 2),
                                                   ensembl_gene_id = gene,
                                                   interaction = paste0(interaction_coef, ":gene"),
                                                   formula_type = type,
                                                   logrank_p = signif(logrank_p, 2),
                                                   wald_p = signif(logrank_p, 2),
                                                   lrt_p = signif(logrank_p, 2),
                                                   formula = form)

                             #Do the same for conf intervals
                             conf = rownames_to_column(as.data.frame(x$conf.int), "feature")
                             conf = conf[conf$feature == paste0(interaction_coef, ":", gene), , drop=F]
                             HR.confint.lower <- signif(conf[,"lower .95"], 2)
                             HR.confint.upper <- signif(conf[,"upper .95"],2)
                             int$HR <- paste0(int$HR, " (",
                                              HR.confint.lower, "-", HR.confint.upper, ")")
                             res<- int %>% select(ensembl_gene_id, beta, HR, p.value, everything())
                             res <- res %>% rename(`HR (95% CI for HR)` = HR)
                             res <- res %>% select(-feature, -`se(beta)`, -z) #redundant
                             res <- as.vector(res[1,,drop=T])

                             return(res)
                             })
    #res_df <- t(as.data.frame(gene_results, check.names = FALSE))

    res_df <- bind_rows(gene_results)
    #return(res_df)
    }

  df = res_df %>%
    left_join(., anno_df, by="ensembl_gene_id") %>%
      mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
      select(gene_name, fdr, beta,`HR (95% CI for HR)`, p.value, gene_type, description, everything()) %>%
      arrange(fdr, p.value)

  return(df)
}


## Apply formula per gene

genewise_cox_interactions <- function(gene_list, time, event, data = coxdata, type="full",
                                      interaction_coef = "ppbc_inv", show_runtime=T){
  stopifnot(type == "full") #"interaction_only" Not yet implemented
  stopifnot(time %in% colnames(data))
  stopifnot(event %in% colnames(data))

  start <- Sys.time()

  #time = "months_overall_survival"
  #event = "overall_survival"
  #data = coxdata
  study_levels = levels(data$study_group)

  clin_cov = 'age_at_diagnosis+year_of_diagnosis+stage+grade+surgery+radiotherapy+hormonetherapy+chemotherapy+herceptin+PAM50+'

  formula = paste0('Surv(time=',time,', event=',event,')~',clin_cov)

  if(type == "full"){
    formula = paste0(formula,paste(study_levels, collapse="+"),"+")
    if(show_runtime){print(formula)}
    gene_formulas <-  sapply(gene_list, function(x) as.formula(paste(formula, paste0(levels(coxdata$study_group),"*", x, collapse = "+"))))

    gene_models <- lapply(gene_formulas, function(x){coxph(x, data = data)})

    res <- extract_interaction_cox_results(gene_models, type=type, interaction_coef = interaction_coef,
                                        time = time, event = event) #Comment out to test
  }

  end <- Sys.time()

  if(show_runtime){print(end-start)}
  #return(res)

  return(res)
}


#Results directory
resDir = here("data", "Rds")
dir.create(resDir, showWarnings = F)
stopifnot(file.exists(resDir))


#### OS ####

print("Genes start at position 25")
print(colnames(coxdata)[25])

resPath = file.path(resDir, "12b_interaction_full_os.Rds")

if(file.exists(resPath) == F| overwrite == T){
  print("Interaction ppbc_inv:gene for OS, full formula")
  res <- genewise_cox_interactions(gene_list = colnames(coxdata)[25:ncol(coxdata)],
                                   time = "months_overall_survival",
                                   event = "overall_survival",
                                   data = coxdata,
                                   type="full",
                                   interaction_coef = "ppbc_inv",
                                   show_runtime=T)
  saveRDS(res, resPath)
}


### DRS ####

resPath = file.path(resDir, "12b_interaction_full_drs.Rds")

if(file.exists(resPath) == F| overwrite == T){
  print("Interaction ppbc_inv:gene for DRS, full formula")
  res <- genewise_cox_interactions(gene_list = colnames(coxdata)[25:ncol(coxdata)],
                                   time = "months_to_drs",
                                   event = "distant_recurrence",
                                   data = coxdata,
                                   type="full",
                                   interaction_coef = "ppbc_inv",
                                   show_runtime=T)
  saveRDS(res, resPath)
}

