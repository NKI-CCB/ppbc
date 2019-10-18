library(DESeq2)
library(apeglm)
library(here)
library(ggrepel)
library(ComplexHeatmap)
library(openxlsx)
library(scrime)
library(RColorBrewer)
library(circlize)
library(tidyverse)

#Script version of notebook 7, without tumor purity

# Load and pre-process data


dds <- readRDS(file = here("data/Rds/05_dds_PAM50_batch.Rds"))

print(design(dds))

gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>% distinct()



## Filtering
#Minimum threshold is a nonzero count in at least a 3rd of all samples

keep <- rowSums(counts(dds)!=0) >= ceiling(ncol(dds)/3)
table(keep)


#Keep the rest
dds <- dds[keep,]



## Retrieve immune genes

#For downstream analyses, we'd like to know which of the differentially expressed genes are immunologically relevant.

#An immune gene is defined as follows:

#Either immune/immuno/interleukin is part of the gene name and description OR the gene is part of the the ImmPort database. See https://www.innatedb.com/redirect.do?go=resourcesGeneLists

#The Immunology Database and Analysis Portal (ImmPort) system was developed under the Bioinformatics Integration Support Contract (BISC) Phase II by the Northrop Grumman Information Technology Health Solutions team for the NIH, NIAID, and DAIT. The principal investigator of the BISC project is Dr. Richard Scheuermann at University of Texas Southwestern Medical Center. The list of immunologically related genes in ImmPort is a collection of ~6,000 human genes, which was formed with the goal of retrieving all genes that have immune system-related functions. This list was generated using automatic searches of EntrezGene and Gene Ontology records using immunology-related keywords. The list was then manually curated by immunology experts examining various literature sources.

immune_gene_list <- read_csv(here("data", "external","gene-sets","InnateDB_genes.csv"))
head(immune_gene_list)


## Helper functions

source(here("src", "deseq_report_functions.R"))


## Variance stabilizing transformation

#We want a transformation that is fully blind to the experimental design, but still uses the faster sampling method in the vst wrapper for varianceStabilizingTransformation.

#This homoskedastic dataset will be used later on for heatmap visualization.
#Previously performed in notebook 6.


vsd = readRDS(here("data","Rds","06_vsd_standardfilter.Rds"))



# Diffex: non-prbc as reference group

#In this call to DESeq, we're focusing in group-by-group contrasts and using a Wald test instead of a likelihood ratio test.

#Ensure that non-prbc is the reference level

dds$study_group = relevel(dds$study_group, ref="non_prbc")

print(levels(dds$study_group))


#Perform differential expression analysis.


start = Sys.time()
dds <- DESeq(dds)
end = Sys.time()
print(end-start)

saveRDS(dds,here("data/Rds/07_dds_apeglm_nonprbcref_standard.Rds"))

# Diffex: prbc as reference level

#From the manual:

#Although apeglm cannot be used with contrast, we note that many designs can be easily rearranged such that what was a contrast becomes its own coefficient. In this case, the dispersion does not have to be estimated again, as the designs are equivalent, up to the meaning of the coefficients. Instead, one need only run nbinomWaldTest to re-estimate MLE coefficients – these are necessary for apeglm – and then run lfcShrink specifying the coefficient of interest in resultsNames(dds). We give some examples below of producing equivalent designs for use with coef. We show how the coefficients change with model.matrix, but the user would, for example, either change the levels of dds$condition or replace the design using design(dds)<-, then run nbinomWaldTest followed by lfcShrink.


dds$study_group = relevel(dds$study_group, ref="prbc")
print(levels(dds$study_group))

start = Sys.time()
dds = nbinomWaldTest(dds)
end = Sys.time()
print(end - start)

saveRDS(dds,here("data/Rds/07_dds_apeglm_prbcref_standard.Rds"))


# Diffex: Lac as reference level

#We have only one comparison remaining: involution vs lac. For this we need to set lactation as the reference group.

dds$study_group = relevel(dds$study_group, ref="ppbc_lac")
print(levels(dds$study_group))


start = Sys.time()
dds = nbinomWaldTest(dds)
end = Sys.time()
print(end - start)

saveRDS(dds,here("data/Rds/07_dds_apeglm_lacref_standard.Rds"))
