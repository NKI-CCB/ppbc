
#export OMP_NUM_THREADS=1 #Run in terminal prior to invoking script, then run script directly from terminal
#Sys.setenv(OMP_NUM_THREADS=1) #This doesn't work from within the script

#This is the R script version of notebook 6 to be run in the background

#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
#install_github("mikelove/DESeq2")
#install_github("azhu513/apeglm")

library(DESeq2)
library(apeglm)
library(here)
library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(openxlsx)
library(scrime)
library(RColorBrewer)

# Load and pre-process data

sink("testingparalellel.txt")

dds <- readRDS(file = here("data/Rds/05_dds_PAM50_batch.Rds"))

print(design(dds))

gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>% distinct()



## Filtering

#Require non-zero counts in 1/3 of the dataset.


keep <- rowSums(counts(dds)!=0) >= ceiling(ncol(dds)/3)
table(keep)
dds <- dds[keep,]



## Retrieve immune genes

#For downstream analyses, we'd like to know which of the differentially expressed genes are immunologically relevant.
#An immune gene is defined as follows:
#Either immune/immuno/interleukin is part of the gene name and description OR the gene is part of the the [ImmPort database](https://www.innatedb.com/redirect.do?go=resourcesGeneLists).
#The Immunology Database and Analysis Portal (ImmPort) system was developed under the Bioinformatics Integration Support Contract (BISC) Phase II by the Northrop Grumman Information Technology Health Solutions team for the NIH, NIAID, and DAIT. The list of immunologically related genes in ImmPort is a collection of ~6,000 human genes, which was formed with the goal of retrieving all genes that have immune system-related functions. This list was generated using automatic searches of EntrezGene and Gene Ontology records using immunology-related keywords. The list was then manually curated by immunology experts examining various literature sources.

immune_gene_list <- read_csv(here("data", "external","gene-sets","InnateDB_genes.csv"))
head(immune_gene_list)

## Helper functions

source(here("src", "deseq_report_functions.R"))

## Variance stabilizing transformation

#We want a transformation that is fully blind to the experimental design, but still uses the faster sampling method in the vst wrapper for varianceStabilizingTransformation.
#This homoskedastic dataset will be used later on for heatmap visualization.

#vsd = blind_vst(dds)
#saveRDS(vsd, here("data","Rds","06_vsd_standardfilter.Rds"))
vsd = readRDS(here("data","Rds","06_vsd_standardfilter.Rds"))

# Likelihood ratio test

#The likelihood ratio test (LRT) is an ANOVA-like test design to check for genes which are significantly differentially expressed in at least one group. From the manual: "The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model, in which some of the terms of the full model are removed. The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero."
#Our design formula is as follows (tumor_purity has been removed on the basis of notebook 5):

print(design(dds))

#We will test it versus a formula that lacks the study group for genes of interest.

start <- Sys.time()
ddsLRT = DESeq(dds, test="LRT", reduced= ~ batch + PAM50)
end <- Sys.time()
print(end-start)

saveRDS(ddsLRT, here("data/Rds/06_dds_testingparallel.Rds"))

sink()