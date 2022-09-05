library(DESeq2)
library(edgeR)
library(tidyverse)
library(here)

resDir = here("data/rnaseq/processed")

#### Load data ----

#Metadata containing clinical covariates, as defined in notebook 4
survdata <- readRDS(here("data/rnaseq/interim/04_survdata.Rds"))
#Exclude unnecessary columns, see notebook 4 for details
survdata <- survdata %>% select(-reason_death, -year_diagnosis)
print(colnames(survdata))

#Gene annotation
gx_annot <- read_tsv(here("data/rnaseq/metadata/01_gene_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id,
                               gene_name, gene_type,
                               description = gene_description) %>% distinct()

#Gene expression
dds = readRDS(here("data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds"))

#### Prepare cox data ----

#Extract raw counts
genEx <- assay(dds)

#Subset for only those samples with enough survival data
print(paste("Samples in gene expression matrix gene expression matrix:", ncol(genEx)))
genEx <- genEx[,colnames(genEx) %in% survdata$sample_name]
print(paste("Samples with adequate metadata:", ncol(genEx)))

#TMM and log transform.
#Anecdotally, TMM does a better job than VST at dealing with outliers

dge <- DGEList(genEx, samples = survdata, group = survdata$study_group)

normTMMlog2 <- function(object){
  object = calcNormFactors(object, method="TMM")
  object = cpm(object, log=T, normalized.lib.sizes=T)
  return(object)
}

mat <- normTMMlog2(dge)

#Combine into single dataframe for Cox
stopifnot(identical(rownames(t(mat)), survdata$sample_name))
coxdata = cbind(survdata, t(mat))

head(coxdata[,1:30])

saveRDS(coxdata, here(resDir, "12_coxdata.Rds"))

#### Involution samples only ----

#Subset for only those samples with enough survival data
print("All samples")
print(ncol(genEx))

print("Involution samples:")
invEx <- genEx[,colnames(genEx)[str_detect(colnames(genEx), "inv")]]
ncol(invEx)

invmeta <- survdata %>% filter(sample_name %in% colnames(invEx))

#TMM and log transform

invdge <- DGEList(invEx, samples = invmeta, group = invmeta$study_group)
invmat <- normTMMlog2(invdge)

#Combine into single dataframe for Cox
stopifnot(identical(rownames(t(invmat)), invmeta$sample_name))

invdata = cbind(invmeta, t(invmat))
inv_gene_col = which(colnames(invdata)=="ENSG00000000003")

print(paste("Gene columns from", inv_gene_col, "to", ncol(invdata), "of involution-only coxdata"))
print(head(invdata[,1:30]))
saveRDS(invdata, file.path(resDir, "12_invdata.Rds"))
