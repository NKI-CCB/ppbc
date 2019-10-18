library(DESeq2)
library(apeglm)
library(here)
library(tidyverse)
library(ggrepel)
library(ComplexHeatmap)
library(openxlsx)
library(scrime)
library(RColorBrewer)

#As notebook 8, but without tumor purity

# Load and pre-process data


dds <- readRDS(file = here("data/Rds/05_dds_PAM50_batch.Rds"))

gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>% distinct()



## Filtering

#Minimum threshold is a nonzero count in at least a 3rd of all samples


keep <- rowSums(counts(dds)!=0) >= ceiling(ncol(dds)/3)
table(keep)

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


## Heatmap colors

#Created in notebook 6

all_colors = readRDS(here("data","Rds", "06_heatmap_colors.Rds"))
study_colors = all_colors$study_colors
pam_colors = all_colors$pam_colors
gene_colors = all_colors$gene_colors


# One vs rest formula levels
print(design(dds))

dds$inv_vs_rest  = colData(dds) %>% as.data.frame() %>% mutate(inv_vs_rest = if_else(study_group == "ppbc_inv", "ppbc_inv", "rest")) %>%
  pull(inv_vs_rest) %>% factor(levels=c("rest", "ppbc_inv"))

dds$prbc_vs_rest  = colData(dds) %>% as.data.frame() %>% mutate(prbc_vs_rest = if_else(study_group == "prbc", "prbc", "rest")) %>%
  pull(prbc_vs_rest) %>% factor(levels=c("rest", "prbc"))

dds$lac_vs_rest  = colData(dds) %>% as.data.frame() %>% mutate(lac_vs_rest = if_else(study_group == "ppbc_lac", "ppbc_lac", "rest")) %>%
  pull(lac_vs_rest) %>% factor(levels=c("rest", "ppbc_lac"))

dds$nonprbc_vs_rest  = colData(dds) %>% as.data.frame() %>% mutate(nonprbc_vs_rest = if_else(study_group == "non_prbc", "non_prbc", "rest")) %>%
  pull(nonprbc_vs_rest) %>% factor(levels=c("rest", "non_prbc"))

table(dds$inv_vs_rest, dds$study_group) %>% print()
table(dds$prbc_vs_rest, dds$study_group) %>% print()
table(dds$lac_vs_rest, dds$study_group) %>% print()
table(dds$nonprbc_vs_rest, dds$study_group) %>% print()



# Diffex: involution vs rest

#In this call to DESeq, we're focusing one one vs rest contrasts and using a Wald test instead of a likelihood ratio test.


dds$group = dds$inv_vs_rest

print(levels(dds$group))

design(dds) = ~batch + PAM50 + group

print(design(dds))

#Perform differential expression analysis.

start = Sys.time()
dds <- DESeq(dds)
end = Sys.time()
print(end-start)

saveRDS(dds,here("data/Rds/08_dds_inv_vs_rest_standard.Rds"))

# Diffex: Pregnancy vs rest

dds$group = dds$prbc_vs_rest
levels(dds$group) %>% print()

#Perform differential expression analysis.

#From the manual:
#Although apeglm cannot be used with contrast, we note that many designs can be easily rearranged such that what was a contrast becomes its own coefficient.
#In this case, the dispersion does not have to be estimated again, as the designs are equivalent, up to the meaning of the coefficients.
#Instead, one need only run nbinomWaldTest to re-estimate MLE coefficients – these are necessary for apeglm – and then run lfcShrink specifying the coefficient of interest in resultsNames(dds). We give some examples below of producing equivalent designs for use with coef. We show how the coefficients change with model.matrix, but the user would, for example, either change the levels of dds$condition or replace the design using design(dds)<-, then run nbinomWaldTest followed by lfcShrink.


start = Sys.time()
dds <- nbinomWaldTest(dds)
end = Sys.time()
print(end-start)

saveRDS(dds,here("data/Rds/08_dds_prbc_vs_rest_standard.Rds"))

# Diffex: Lac vs rest

dds$group = dds$lac_vs_rest
levels(dds$group) %>% print()

#Perform differential expression analysis.

start = Sys.time()
dds <- nbinomWaldTest(dds)
end = Sys.time()
print(end-start)

saveRDS(dds,here("data/Rds/08_dds_lac_vs_rest_standard.Rds"))

# Diffex: Non-prbc vs rest

dds$group = dds$nonprbc_vs_rest
print(levels(dds$group))


#Perform differential expression analysis.


start = Sys.time()
dds <- nbinomWaldTest(dds)
end = Sys.time()
print(end-start)

saveRDS(dds,here("data/Rds/08_dds_nonprbc_vs_rest_standard.Rds"))
