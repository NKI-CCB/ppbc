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
resDir = here::here("data", "rnaseq", "interim")
dir.create(resDir, showWarnings = F)
stopifnot(file.exists(resDir))

shrinkRes <- function(dds, contrast, type="apeglm"){

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
dds <- readRDS(here("data/rnaseq/interim/05b_dds_filtered.Rds"))
print(paste(nrow(dds), "genes in dataset"))

design(dds) <- ~batch + PAM50 + study_group

print("Design formula:")
print(design(dds))

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
#And that involution is always the first group, where possible
for (i in 1:nrow(combos)){
  first <- combos[i,1]
  second <- combos[i,2]
  if (first == "npbc" | second == "ppbcpw"){
    combos[i,1] <- second
    combos[i, 2] <- first
  }
}
combos <- as.data.frame(combos)
colnames(combos) <- c("primary", "reference")
print("Pairwise comparisons:")
print(combos)

#Apeglm files:
#combos %>%
#  mutate(comp = paste(primary, reference, sep = "_vs_")) %>%
#  pull(comp)

#Keep unanalyzed copy of dds separate
ddsR <- dds

#Pairwise loop
for (ref in unique(combos$reference)){
  
  pw_dds_file <- file.path(resDir, paste0("07_dds_pairwise_ref_", ref, ".Rds"))
  
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
    
    pw_ape_file <- file.path(resDir, paste0("07_ape_", primary, "_vs_", ref, ".Rds"))
    
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
