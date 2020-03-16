# Differential expression jobs to be displayed in notebook 17.

rm(list = ls())

library(here)
library(tictoc)
library(DESeq2)
library(apeglm)
library(tidyverse)

#### Important ----

#To speed up DESeq
#Open the bash terminal and log into harris/darwin
#Run this command:
#export OMP_NUM_THREADS=1
#Then run the script from the console

#### Setup ----

#Whether to overwrite existing results
overwrite <- T

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
dds <- readRDS(here("data","Rds","17_dds_receptorstatus.Rds"))
print(paste(nrow(dds), "genes in dataset"))

design(dds) <- ~ER + HER2 + study_group

print("Design formula:")
print(design(dds))

### LRT ----

lrtPath = file.path(resDir, "17_dds_LRT.Rds")

if(file.exists(lrtPath) == F| overwrite == T){
  print("Likelihood ratio test vs reduced formula ~ER + HER2")
  tic("LRT")
  ddsLRT <- DESeq(dds, test="LRT", reduced= ~ ER + HER2)
  saveRDS(ddsLRT, lrtPath)
  toc()
}

#### Pairwise comparisons ----

#In the first iteration, we call DESeq2()
#Afterwards we call nbinomWaldTest because we are merely relevelling the design 
#and as such, do not need to re-estimate the dispersions

#From the manual:
#Although apeglm cannot be used with contrast, we note that many designs can be easily 
#rearranged such that what was a contrast becomes its own coefficient.
#In this case, the dispersion does not have to be estimated again, as the designs are equivalent, up to the meaning of the coefficients.
#Instead, one need only run nbinomWaldTest to re-estimate MLE coefficients – these are necessary for apeglm – 
#and then run lfcShrink specifying the coefficient of
#interest in resultsNames(dds.nolac). We give some examples below of producing equivalent designs for use with coef.
#We show how the coefficients change with model.matrix, but the user would, for example, 
#either change the levels of dds$condition or replace the design using design(dds)<-,
#then run nbinomWaldTest followed by lfcShrink.

overwrite <- T

#Set up the possible comparisons
combos <- t(combn(levels(dds$study_group), 2))

#Ensure that nulliparous is always the reference (2nd group) where possible
#And that involution is always the first group where possible
for (i in 1:nrow(combos)){
  first <- combos[i,1]
  second <- combos[i,2]
  if (first == "non_prbc" | second == "ppbc_inv"){
    combos[i,1] <- second
    combos[i, 2] <- first
  }
}
combos <- as.data.frame(combos)
colnames(combos) <- c("primary", "reference")
print("Pairwise comparisons:")
print(combos)

#Keep unanalyzed copy of dds separate
ddsR <- dds

#Pairwise loop
for (ref in unique(combos$reference)){
  
  pw_dds_file <- file.path(resDir, paste0("17_dds_pairwise_ref_", ref, ".Rds"))
  
  if (file.exists(pw_dds_file) == F | overwrite == T){
    
    #Relevel for current reference
    print(paste("Pairwise comparisons with ref:", ref))
    ddsR$study_group <- relevel(ddsR$study_group, ref = as.character(ref))
    print(design(ddsR))
    print(paste("Levels:", paste(levels(ddsR$study_group), collapse= ", ")))
    
    print("Results dds at:")
    print(pw_dds_file)
    
    #Differential expression
    #If we haven't run DESeq2 yet, or if this is the first time with this design, run it
    if(ref == unique(combos$reference)[1] | is.null(dispersions(ddsR))){
      tic(paste("DESeq2 with ref", ref))
      ddsR <- DESeq2::DESeq(ddsR)
      saveRDS(object = ddsR, file = pw_dds_file)
      toc()
    } else {
      tic(paste("NbinomWaldTest with ref", ref))
      ddsR <- DESeq2::nbinomWaldTest(ddsR)
      toc()
      saveRDS(object = ddsR, file = pw_dds_file)
    }
    

  }
  
  #Select comparisons that involve the current reference
  comps <- combos[combos$reference==ref,,drop=F]
  
  for (primary in unique(comps$primary)){
   
    pw_ape_file <- file.path(resDir, paste0("17_ape_", primary, "_vs_", ref, ".Rds"))
    
    if (file.exists(pw_ape_file) == F | overwrite == T){
      print("Reloading results from file:")
      print(pw_dds_file)
      ddsR <- readRDS(pw_dds_file)
      print(paste("Apeglm fold shrinkage for", primary, "vs", ref))
      print("Apeglm results at:")
      print(pw_ape_file)
      
      #Fold shrinkage with apeglm
      tic("apeglm")
      ape = shrinkRes(ddsR,  contrast=c("study_group", primary, ref))
      saveRDS(object = ape, file = pw_ape_file)
      toc()
      
    }
  }
}

#### One vs rest comparisons ----

overwrite <- T

ovr_comps <- colnames(colData(dds))[str_detect(colnames(colData(dds)), "vs_rest")]

print("One vs rest comparisons:")
print(paste(ovr_comps, collapse = ", "))

#Keep unanalyzed copy separate
ddsOVR <- dds

for (ovr in ovr_comps){

  design(ddsOVR) <- as.formula(paste("~ ", paste(c("ER", "HER2", ovr), collapse= "+")))
  ovr_dds_file <- file.path(resDir, paste0("17_dds_ovr_", ovr, ".Rds"))

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
  
  ovr_ape_file <- file.path(resDir, paste0("17_ape_", ovr, ".Rds"))
  
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


