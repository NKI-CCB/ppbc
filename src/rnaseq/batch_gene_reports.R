library(here)
library(rmarkdown)
library(tidyverse)

source(here("src/utils/parse_args.R"))

argdf <- retrieve_args()

print(argdf)

#List of genes from which to generate reports
genes_to_report = read_csv(
  here(filter(argdf, argname == "genes")$argval)
)

outdir=here(filter(argdf, argname == "outdir")$argval)
print(paste("Outdir:", outdir))
dir.create(here(outdir), showWarnings = F)

#Allows naming report with both gene symbol and ensembl_id
bx_annot <- readRDS(here("data/rnaseq/processed/bx_annot.Rds"))
dds = readRDS(here("data/rnaseq/processed/08_dds_ovr_inv_vs_rest.Rds"))
genes_passing_filter = rownames(dds)

gene_lookup <- function(gene, dict=bx_annot){

  if(gene %in% dict$gene_name){
    gene_name = gene
    ensembl_gene_id = unique(dict[dict$gene_name == gene_name,]$ensembl_gene_id)
  } else if(gene %in% dict$ensembl_gene_id) {
    ensembl_gene_id = gene
    gene_name = unique(dict[dict$ensembl_gene_id == ensembl_gene_id,]$gene_name)
  } else {
    stop("Provided gene is not a valid ensembl id or gene name")
  }

  uni=unique(dict[dict$ensembl_gene_id == ensembl_gene_id[1],]$uniprot_id)
  ent=unique(dict[dict$ensembl_gene_id == ensembl_gene_id[1],]$entrez_id)

  list(gene_name=gene_name,
       ensembl_gene_id = ensembl_gene_id,
       uniprot_id=uni,
       entrez_id=ent)
}

#Project directory and parameters
projDir <- here()
Rmd <- file.path(projDir, "reports/rnaseq/17_gene_report_template.Rmd")
stopifnot(file.exists(Rmd))
overwrite = F

#Loop over every gene and create report
for(gene in genes_to_report$Gene){

  print(paste("Generating report for", gene))

  #If provided input is ensembl ID, lookup the name of the gene
  #Title will be gene name or gene name and ens ID, if ens ID is provided
  #If ensembl ID is not provided and gene contains more than one,
  #All matching ensembl IDs will be used in the report
  if (gene %in% bx_annot$ensembl_gene_id){
    gn <- gene_lookup(gene)$gene_name
    ens <- gene
    title=paste(gn, ens, sep=":")
  } else {
    title <- gene
  }

  print(title)

  outfile=file.path(outdir, paste0(title, "_report.pdf"))
  print(paste("Writing report to", outfile))

  if (!file.exists(outfile) | overwrite == T){
    rmarkdown::render(input =  Rmd,
                      output_file = outfile,
                      params = list(gene = gene,
                                    show_functions = F,
                                    show_code = T))
  }
}


