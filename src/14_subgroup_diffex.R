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
overwrite <- T

#Results directory
resDir = here::here("data", "Rds", "subgroup_diffex")
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

#Gene expression
dds = readRDS(here("data/Rds/08_dds_ovr_inv_vs_rest.Rds"))

# Filtering
#Minimum threshold is a nonzero count in at least a 3rd of all samples

keep <- rowSums(counts(dds)!=0) >= ceiling(ncol(dds)/3)
table(keep)

#Keep the rest
dds <- dds[keep,]

#Test dds
#dds = head(dds)

#Remove PAM50 from design, as we'll be splitting the dds by PAM50 instead
design(dds) <- ~ batch + study_group

#### Subgroup comparisons ----
print("Subgroups:")
print(table(dds$PAM50, dds$study_group))
print(table(dds$PAM50, dds$inv_vs_rest))
#print(table(dds$ER, dds$inv_vs_rest))

#Too few samples in PAM50:normal to include
#Too few to do pairwise comparisons with lac as well
subgroups = levels(dds$PAM50)[-5]
print("Subgroups:")
print(subgroups)

comps = tibble(
  group = rep("ppbc_inv", 3),
  ref = c("non_prbc", "prbc", "rest"),
  type = c(rep("pairwise", 2), "inv_vs_rest")
)
print("Comparisons:")
comps


for (subgroup in subgroups){
  
  #subgroup = "Basal" #Testing
  print(paste("Starting analysis for subgroup:", subgroup))
  
  dds_sub = dds[, dds$PAM50 == subgroup]
  dds_sub
  
  for (i in 1:nrow(comps)){
    
    #i = 1 #testing
    group = comps$group[i]
    ref = comps$ref[i]
    type = comps$type[i]
    
    #Generate filename based on subgroup comparison
    print(paste0("Diffex ", group, " vs ", ref, ", subgroup: ", subgroup))
    dds_file=file.path(resDir, paste0("14_dds_", subgroup, "_", group, "_vs_", ref, ".Rds"))
    
    #Apply the correct formula based on the type of comparisons
    if(type == "pairwise"){
      dds_sub$study_group = relevel(dds_sub$study_group, ref=ref)
      design(dds_sub) <- ~ batch + study_group
      column = "study_group"
      print(paste("Pairwise design:", paste0(design(dds_sub), collapse="")))
    } else {
      dds_sub$inv_vs_rest <- relevel(dds_sub$inv_vs_rest, ref=ref)
      design(dds_sub) <- ~ batch + inv_vs_rest
      column = "inv_vs_rest"
      print(paste("Inv vs rest design:", paste0(design(dds_sub), collapse="")))
    }
    
    #Report n samples in comparisons
    subdf = as.data.frame(colData(dds_sub)[,c("PAM50", column)])
    print(paste(sum(subdf[,2]==group), subgroup, "sample in", group))
    print("vs")
    print(paste(sum(subdf[,2]==ref), subgroup, "sample in", ref))
    sum(subdf[,2]==ref)
    
    #Run DESeq
    if(!file.exists(dds_file) | overwrite == T){
      tic("DESeq")
      dds_sub <- DESeq2::DESeq(dds_sub)
      toc()
      saveRDS(dds_sub, dds_file)
      print("dds results at:")
      print(dds_file)
    } 
    
    #Shrink fold changes with apeglm
    ape_file <- file.path(resDir, paste0("14_ape_", subgroup, "_", group, "_vs_", ref, ".Rds"))
    
    if (file.exists(ape_file) == F | overwrite == T){
      print("Continuing with results from file:")
      print(dds_file)
      dds_sub <- readRDS(dds_file)
      print(paste("Apeglm fold shrinkage for", group, "vs", ref))
      
      #Fold shrinkage with apeglm
      tic("apeglm")
      ape = shrinkRes(dds_sub,  contrast=c(column, group, ref))
      toc()
      saveRDS(object = ape, file = ape_file)
      print("Apeglm results at:")
      print(ape_file)
      print("")
    }

  }
}
