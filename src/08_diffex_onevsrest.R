library(here)
library(tictoc)
library(DESeq2)
library(apeglm)
library(tidyverse)

rm(list = ls())

#### Important ----

#To speed up DESeq
#Open the bash terminal and log into harris/darwin
#Run this command:
#export OMP_NUM_THREADS=1
#Then run the script from the console

#### Setup ----

#Whether to overwrite existing results
overwrite <- F

#Results directory
resDir = here::here("data", "Rds")
dir.create(resDir, showWarnings = F)
stopifnot(file.exists(resDir))

shrinkRes <- function(dds, contrast, type="apeglm"){
  require(DESeq2)
  require(tidyverse)
  
  cdf = enframe(colnames(coef(dds)), "coef", "name")
  contrast_name = paste(contrast[1], contrast[2], "vs",contrast[3], sep="_")
  
  if(!contrast_name %in% cdf$name){
    stop("Contrast name not found within coefficients. Do you need to relevel your design formula?")
  }
  
  print("Detecting coefficient corresponding to contrast:")
  
  c = cdf[cdf$name==contrast_name,]
  print(paste("Contrast name:", c$name))
  print(paste("Coefficient:", c$coef))
  
  print("Retrieving results from deseq dataset with contrast...")
  res = results(dds, contrast=contrast)
  
  print(paste("Shrinking log2foldchanges via method", type))
  ape <- lfcShrink(dds, coef = c$coef, res = res, type=type)
  
  return(ape)
}

#### Load data ----

#dds with filtering and one vs rest groups pre-defined, see notebook 8
dds <- readRDS(here("data","Rds","05_dds_PAM50_batch.Rds"))
print(paste(nrow(dds), "genes in dataset"))

# Filtering
#Minimum threshold is a nonzero count in at least a 3rd of all samples

keep <- rowSums(counts(dds)!=0) >= ceiling(ncol(dds)/3)
table(keep)

#Keep the rest
dds <- dds[keep,]

#### One vs rest comparisons ----

dds$inv_vs_rest  = colData(dds) %>% as.data.frame() %>%
  mutate(inv_vs_rest = if_else(study_group == "ppbc_inv", "ppbc_inv", "rest")) %>%
  pull(inv_vs_rest) %>%
  factor(levels=c("rest", "ppbc_inv"))

dds$prbc_vs_rest  = colData(dds) %>% as.data.frame() %>%
  mutate(prbc_vs_rest = if_else(study_group == "prbc", "prbc", "rest")) %>%
  pull(prbc_vs_rest) %>% factor(levels=c("rest", "prbc"))

dds$lac_vs_rest  = colData(dds) %>% as.data.frame() %>%
  mutate(lac_vs_rest = if_else(study_group == "ppbc_lac", "ppbc_lac", "rest")) %>%
  pull(lac_vs_rest) %>% factor(levels=c("rest", "ppbc_lac"))

dds$nonprbc_vs_rest  = colData(dds) %>% as.data.frame() %>%
  mutate(nonprbc_vs_rest = if_else(study_group == "non_prbc", "non_prbc", "rest")) %>%
  pull(nonprbc_vs_rest) %>% factor(levels=c("rest", "non_prbc"))

table(dds$inv_vs_rest, dds$study_group)
table(dds$prbc_vs_rest, dds$study_group)
table(dds$lac_vs_rest, dds$study_group)
table(dds$nonprbc_vs_rest, dds$study_group)


ovr_comps <- colnames(colData(dds))[str_detect(colnames(colData(dds)), "vs_rest")]

print("One vs rest comparisons:")
print(paste(ovr_comps, collapse = ", "))

#Keep unanalyzed copy separate
ddsOVR <- dds

for (ovr in ovr_comps){
  
  design(ddsOVR) <- as.formula(paste("~ ", paste(c("batch", "PAM50", ovr), collapse= "+")))
  ovr_dds_file <- file.path(resDir, paste0("08_dds_ovr_", ovr, ".Rds"))
  
  if(file.exists(ovr_dds_file) ==F | overwrite == T){
    
    print(paste("One vs rest comparison with design:"))
    print(design(ddsOVR))
    print("Results file at:")
    print(ovr_dds_file)
    
    #Run DESeq2 if this is the first run with this design, or if dispersions are lacking
    if(ovr == ovr_comps[1] | is.null(dispersions(ddsOVR))){
      tic(paste("DESeq with:", ovr))
      ddsOVR <- DESeq2::DESeq(ddsOVR)
      saveRDS(object = ddsOVR, file = ovr_dds_file)
      toc()
    } else {
      tic(paste("NbinomWaldTest with:", ovr))
      ddsOVR <- DESeq2::nbinomWaldTest(ddsOVR)
      saveRDS(object = ddsOVR, file = ovr_dds_file)
      toc()
    }
    
  }
  
  ovr_ape_file <- file.path(resDir, paste0("08_ape_ovr_", ovr, ".Rds"))
  
  if(file.exists(ovr_ape_file) == F | overwrite == T){
    
    #Reload the dds to ensure we have the right one
    print("Reloading results from file:")
    print(ovr_dds_file)
    ddsOVR <- readRDS(file = ovr_dds_file)
    
    primary <- levels(colData(ddsOVR)[,ovr])[2]
    print(paste("Apeglm fold shrinkage for", primary, "vs rest"))
    print("Apeglm results at:")
    print(ovr_ape_file)
    
    #Fold shrinkage with apeglm
    tic("apeglm")
    ape = shrinkRes(ddsOVR,  contrast=c(ovr, primary, "rest"))
    saveRDS(object = ape, file = ovr_ape_file)
    toc()
    
  }
  
}