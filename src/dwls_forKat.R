## install dwls algorithm ##
library(devtools)
library(here)

#Installation failed on 3.4.4 and 3.6.1 for reasons too tedious to enumerate
#Eventually working with R 3.5.1
#On darwin: conda create -n R-3.5.1
#devtools::install_github("hhoeflin/hdf5r") #Must use latest version or Seurat cannot be installed
#install.packages('Seurat') #Must be working before dwls can be installed
#devtools::install_bitbucket("yuanlab/dwls", ref="default")
library(DWLS)

source(here("src", "dwls_functions.R"))

target <- read.csv(here("data", "ppbc_samplesxgenes.csv"))
rownames(target) <- target$sample_name
target <- target[,-1]
target <- t(target)
target[1:10,1:10]

refProfiles <- read.table("data/dwls/20200106_refProfiles_STARres0.1CellTypes.csv")
markerGenes <- read.table("data/dwls/20200106_markerGenes_STARres0.1CellTypes.csv")
markerGenes <- markerGenes[,1]
identical(refProfiles, markerGenes)


## target is the matrix with the real tumor samples ##
## refProfiles is the signature matrix ##
## markerGenes is the vector containing all the marker genes ##

## normalize reference profiles and target dataset ##
rt <- as.data.frame(apply(target,2,function(x) x/sum(x))[markerGenes,])
sig <- as.matrix(apply(refProfiles,2,function(x) x/sum(x))[markerGenes,])

## run paralellized dwls on nCores number of cores ##
nCores <- 10
cl <- parallel::makeCluster(nCores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = ncol(rt), style = 3)
progressBar <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progressBar)
results <- foreach(x = c(1:ncol(rt)), .combine=rbind, .options.snow=opts) %dopar% {
  library(DWLS)
  solveDampenedWLS(sig,rt[,x])
}
stopCluster(cl)
rownames(results) <- colnames(rt)
results <- as.data.frame(results)

## calculate some statistics on the deconvolution results ##
results$pearsonCor <- sapply(c(1:nrow(results)),function(x) cor(sig %*% as.numeric(results[x,c(1:ncol(sig))]),rt[,x]))
results$pearsonCorPval <- sapply(c(1:nrow(results)),function(x) cor.test(sig %*% as.numeric(results[x,c(1:ncol(sig))]),rt[,x])$p.value)
results$spearmanCor <- sapply(c(1:nrow(results)),function(x) cor(sig %*% as.numeric(results[x,c(1:ncol(sig))]),rt[,x],method = "spearman"))
results$spearmanCorPval <- sapply(c(1:nrow(results)),function(x) cor.test(sig %*% as.numeric(results[x,c(1:ncol(sig))]),rt[,x],method = "spearman")$p.value)
results$rmse <- sapply(c(1:nrow(results)),function(x) sqrt(mean((sig %*% as.numeric(results[x,c(1:ncol(sig))])-rt[,x])**2)))

## calculate p-value for deconvolution results, based on random tumor samples (same procedure as used in cibersort) ##

## number of permutations ##
perm <- 1000
nCores <- 10
cl <- makeCluster(nCores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = perm, style = 3)
progressBar <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progressBar)
corDist <- foreach ( i = c(1:perm), .combine=c , .options.snow=opts) %dopar% {
  library(DWLS)
  set.seed(i)
  rtlist <- as.list(data.matrix(rt))
  randomMixture <- as.numeric(rtlist[sample(length(rtlist),nrow(sig))])
  randomMixture <- randomMixture/sum(randomMixture)
  cf <- solveDampenedWLS(sig,randomMixture)
  return(cor(sig %*% cf,randomMixture))
}

results$dwlsPval <- sapply(c(1:nrow(results)),function(x) sum(corDist > results$pearsonCor[x]))

## save results ##

write.table(results,"dwls_results.csv")

