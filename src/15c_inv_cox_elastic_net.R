library(here)
library(survival)
library(glmnet)
library(tidyverse)

#As script 15b, but using all genes as input
# Load data 
#A samples x features matrix, where features contain both clinical covariates and log2-transformed, TMM-normalized gene expression data.
coxdata = readRDS(here("data/Rds/12_invdata.Rds"))
#head(coxdata[,1:30])
genecol = which(colnames(coxdata)=="ENSG00000000003")
metacol = 4


#First `r metacol` columns are metadata:
head(coxdata[,1:metacol])

#Clinical covariates and survival/metastasis endpoints are columns 5 through `r genecol-1`.
head(coxdata[,(metacol + 1):(genecol-1)])

#Create separate data frames for metadata and covariates.
metadata <- coxdata %>%
  select(sample_name:PPBC,
         months_of_followup,
         overall_survival,
         months_overall_survival,
         metastasis_at_diagnosis,
         distant_recurrence,
         months_to_drs)
colnames(metadata)


covdf <- coxdata[,1:(genecol-1)]
covdf <- covdf[,!colnames(covdf) %in% colnames(metadata)]
colnames(covdf)
covdf <- data.matrix(covdf)

#Columns `r genecol-1` through `r ncol(coxdata)` are gene expression.

geneEx = coxdata[,genecol:ncol(coxdata)]
head(geneEx[,1:4])

## Gene sets and annotation
gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id, gene_name,
                               gene_type, description = gene_description) %>% distinct()

# Combine only those necessary covariate columns with gene expression
survdata <- cbind(covdf, geneEx)
colnames(survdata)[1:20]
# OS/DRS Time Points

#For family="cox", y should be a two-column matrix with columns named 'time' and 'status'. 
#The latter is a binary variable, with '1' indicating death, and '0' indicating right censored.  

#Overall survival:
time_os = data.frame(time = metadata$months_overall_survival,
                     event = metadata$overall_survival)
rownames(time_os) = metadata$sample_name

#os = as.matrix(os)
time_os = Surv(time_os$time, time_os$event)
head(time_os)

#Metastasis
#When looking at drs, time point 0 (metastatic at diagnosis) is not allowed
met_at_diag <- metadata %>%
  filter(metastasis_at_diagnosis == 1) %>%
  pull(sample_name)

met_at_diag


time_drs = data.frame(time = coxdata$months_to_drs, event = coxdata$distant_recurrence)
rownames(time_drs) = coxdata$sample_name
time_drs = time_drs[!rownames(time_drs) %in% met_at_diag,]

time_drs = time_drs[time_drs$time!=0,]
#drs = as.matrix(drs)
time_drs = Surv(time_drs$time, time_drs$event)
head(time_drs)

survdata.drs = survdata[!(rownames(survdata) %in% met_at_diag),]


# Run Cox elastic net models {.tabset}

#A pentalized cox model to determine which groups of covariates are most predictive of survival.
#Using nested crossvalidation.
nested.cox <- function(x,y,nfolds=10,nfolds_inner=10,
                       n_preds = 5,
                       print_time = T,
                       s="lambda.min",
                       alpha=0.5){
  #Based on nested.cv from the TANDEM package
  start <- Sys.time()
  samples = rownames(x)
  n = length(y)
  foldid = ceiling(sample(1:n)/n * nfolds)
  
  fitted_relative_risk = rep(NA, length(y))
  risklist <- list()
  features = list()
  pred = list()
  fits = list()
  
  for (n in 1:n_preds){
    
    for (i in 1:nfolds) {
          
        ind = foldid == i
        x_train = x[!ind, ]
        y_train = y[!ind]
        x_test = x[ind, ]
        n_i = length(y_train)
        foldid_i = ceiling(sample(1:n_i)/n_i * nfolds_inner)
        fit = glmnet::cv.glmnet(x_train, y_train, family = "cox", alpha = alpha,
                                foldid = foldid_i)
        fits[[length(fits) + 1]] = fit
      
        features[[length(features) + 1]] = coef(fit, s=s)[,1] %>% enframe() %>% filter(value!=0 & name != "(Intercept)")
        fitted_relative_risk[ind] = glmnet::predict.cv.glmnet(fit, 
                newx = x_test, s = s, alpha=alpha, type="response")

    }
    
    relrisk <- tibble(samples = samples, fitted_relative_risk = fitted_relative_risk)
    risklist[[length(risklist)+1]] <- relrisk
    #pred[[length(pred)+1]] <- list(features = features, relrisk = relrisk, fits = fits)
    pred$features <- features
    pred$relrisk <- risklist
    pred$fits <- fits
    print(n)

  }

  end <- Sys.time()
  if (print_time){print(end-start)}
  return(pred)
}

#Ensure matrix class
survdata <- as.matrix(survdata)
survdata.drs <- as.matrix(survdata.drs)

## Overall survival - all genes as input
resDir <- here("data/Rds")
stopifnot(dir.exists(resDir))
set.seed(123)
glm.os <- nested.cox(x = survdata, y = time_os, s = "lambda.min", alpha = 0.5)
saveRDS(glm.os, file.path(resDir,"15c_inv_glm_os_all.Rds"))

## Distant recurrence - all genes as input
set.seed(123)
glm.drs <- nested.cox(x = survdata.drs, y = time_drs, s = "lambda.min", alpha = 0.5)
saveRDS(glm.drs, file.path(resDir,"15c_inv_glm_drs_all.Rds"))

