library(DESeq2)
library(limma)
library(edgeR)
library(tidyverse)
library(doParallel)
library(flexgsea) #https://github.com/NKI-CCB/flexgsea-r
library(here)

dir.create(here("results", "flexgsea", "logs"), showWarnings = F, recursive = T)
sink(here("results", "flexgsea", "logs", "flexgsea_limma_log.txt"))

#Set number of permutations for flexgsea
nperm = 1000 #5 for testing, 1000 for real thing
print(paste("Number of permutations:", nperm))


# Main comparisons of interest
# Involution vs nulliparous
# Pregnant vs nulliparous
# Involution vs all groups (rest)

registerDoParallel(cores=6)

####Load data----

dds = readRDS(here("data", "Rds", "05_dds_PAM50_tumorpurity_batch.Rds"))
gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>% distinct()

#Load gene signatures

gene_signatures = list(
  hallmark = read_gmt(file = here("data","external","gmt","h.all.v6.2.symbols.gmt")), # 50 gene overview
  go_c5 = read_gmt(file = here("data","external","gmt","c5.all.v6.2.symbols.gmt")), # Gene ontology
  oncogene_c6 = read_gmt(file = here("data","external","gmt","c6.all.v6.2.symbols.gmt")), #Oncogenic gene signatures
  immune_c7 = read_gmt(file = here("data","external","gmt","c7.all.v6.2.symbols.gmt")), #Immunological signatures
  canonpath_c2 = read_gmt(file = here("data","external","gmt","c2.cp.v6.2.symbols.gmt")), #Canonical pathways
  motif_tf_mir_c3 = read_gmt(file = here("data","external", "gmt","c3.all.v6.2.symbols.gmt")), #Targets of TFs and miRs
  comp_canc_c4 = read_gmt(file = here("data","external","gmt","c4.all.v6.2.symbols.gmt")) #Computational gene sets from cancer gene neighborhoods and cancer modules
)

signature_files = enframe(c(
  hallmark = here("data","external","gmt","h.all.v6.2.symbols.gmt"), # 50 gene overview
  go_c5 = here("data","external","gmt","c5.all.v6.2.symbols.gmt"), # Gene ontology
  oncogene_c6 = here("data","external","gmt","c6.all.v6.2.symbols.gmt"), #Oncogenic gene signatures
  immune_c7 =  here("data","external","gmt","c7.all.v6.2.symbols.gmt"), #Immunological signatures
  canonpath_c2 =  here("data","external","gmt","c2.cp.v6.2.symbols.gmt"), #Canonical pathways
  motif_tf_mir_c3 = here("data","external", "gmt","c3.all.v6.2.symbols.gmt"), #Targets of TFs and miRs
  comp_canc_c4 = here("data","external","gmt","c4.all.v6.2.symbols.gmt")), #Computational gene sets from cancer gene neighborhoods and cancer modules
  "signature_name", "file")

####Filtering ----

#For limma, higher count filtering is crucial, lest the variance be underrestimated when the mean is low
#The cpm must be at least 1 in 1/3 of the dataset
keep <- rowSums(fpm(dds,robust=F) > 1 ) >= ceiling(ncol(dds)/3)
print("Using limma-style low count filter")
print("Genes with a cpm of at least 1 in 1/3 of the dataset:")
print(table(keep))
dds <- dds[keep,]

#### Define comparisons ----

# Add one vs rest comparison levels
dds$inv_vs_rest  = colData(dds) %>% as.data.frame() %>%
  mutate(inv_vs_rest = if_else(study_group == "ppbc_inv", "inv", "rest")) %>%
  pull(inv_vs_rest) %>% factor(levels=c("rest", "inv"))

dds$prbc_vs_rest  = colData(dds) %>% as.data.frame() %>%
  mutate(prbc_vs_rest = if_else(study_group == "prbc", "prbc", "rest")) %>%
  pull(prbc_vs_rest) %>% factor(levels=c("rest", "prbc"))

dds$lac_vs_rest  = colData(dds) %>% as.data.frame() %>%
  mutate(lac_vs_rest = if_else(study_group == "ppbc_lac", "lac", "rest")) %>%
  pull(lac_vs_rest) %>% factor(levels=c("rest", "lac"))

dds$nonprbc_vs_rest  = colData(dds) %>% as.data.frame() %>%
  mutate(nonprbc_vs_rest = if_else(study_group == "non_prbc", "nonprbc", "rest")) %>%
  pull(nonprbc_vs_rest) %>% factor(levels=c("rest", "nonprbc"))

print("New one vs rest comparisons")
print(table(dds$inv_vs_rest, dds$study_group))
print(table(dds$prbc_vs_rest, dds$study_group))
print(table(dds$lac_vs_rest, dds$study_group))
print(table(dds$nonprbc_vs_rest, dds$study_group))

# Expand study group into pairwise comparison columns
# Due to the way flexgsea handles the coeficients, the selection is simplest when the comparison of interest is a factor with two levels
# The "other" level will be removed downstream
dds$inv_vs_nonprbc  = colData(dds) %>% as.data.frame() %>%
  mutate(inv_vs_nonprbc = case_when(
    study_group == "ppbc_inv" ~ "inv",
    study_group == "non_prbc" ~ "nonprbc",
    TRUE ~ "other")) %>%
  pull(inv_vs_nonprbc) %>% factor(levels=c("nonprbc", "other", "inv"))

dds$prbc_vs_nonprbc = colData(dds) %>% as.data.frame() %>%
  mutate(prbc_vs_nonprbc = case_when(
    study_group == "prbc" ~ "prbc",
    study_group == "non_prbc" ~ "nonprbc",
    TRUE ~ "other")) %>%
  pull(prbc_vs_nonprbc) %>% factor(levels=c("nonprbc", "other", "prbc"))

dds$inv_vs_prbc = colData(dds) %>% as.data.frame() %>%
  mutate(inv_vs_prbc = case_when(
    study_group == "ppbc_inv" ~ "inv",
    study_group == "prbc" ~ "prbc",
    TRUE ~ "other")) %>%
  pull(inv_vs_prbc) %>% factor(levels=c("prbc", "other", "inv"))

dds$inv_vs_lac = colData(dds) %>% as.data.frame() %>%
  mutate(inv_vs_lac = case_when(
    study_group == "ppbc_inv" ~ "inv",
    study_group == "ppbc_lac" ~ "lac",
    TRUE ~ "other")) %>%
  pull(inv_vs_lac) %>% factor(levels=c("lac", "other", "inv"))

dds$lac_vs_nonprbc = colData(dds) %>% as.data.frame() %>%
  mutate(lac_vs_nonprbc = case_when(
    study_group == "ppbc_lac" ~ "lac",
    study_group == "non_prbc" ~ "nonprbc",
    TRUE ~ "other")) %>%
  pull(lac_vs_nonprbc) %>% factor(levels=c("nonprbc", "other", "lac"))

dds$lac_vs_prbc = colData(dds) %>% as.data.frame() %>%
  mutate(lac_vs_prbc = case_when(
    study_group == "ppbc_lac" ~ "lac",
    study_group == "prbc" ~ "prbc",
    TRUE ~ "other")) %>%
  pull(lac_vs_prbc) %>% factor(levels=c("prbc", "other", "lac"))

print("New pairwise comparisons (other will be removed downstream)")
print(table(dds$inv_vs_nonprbc, dds$study_group))
print(table(dds$prbc_vs_nonprbc, dds$study_group))
print(table(dds$inv_vs_prbc, dds$study_group))
print(table(dds$inv_vs_lac, dds$study_group))
print(table(dds$lac_vs_nonprbc, dds$study_group))
print(table(dds$lac_vs_prbc, dds$study_group))

### Convert to gene symbol count matrix ----

#Extract raw counts
counts <- counts(dds, normalized = F )

#Gmt files never have Ensembl gene ids, so we'll need to convert to gene symbols.

df <- right_join(gx_annot[,c("ensembl_gene_id","gene_name")],rownames_to_column(as.data.frame(counts),"ensembl_gene_id"),
                 by="ensembl_gene_id")

print(paste("Number of ensembl IDs with no corresponding gene symbol:", sum(is.na(df$gene_name))))
df <- df[!is.na(df$gene_name),]

#When dealing with duplicated gene names, take the sum expression

summarize_expression_duplicate_ids <- function(mat, id_column, f=colMeans, final_gene_symbol_colname="GeneSymbol"){
  require(dplyr)

  #Easiest way to write functions with dplyr is to standardize the column name

  input = mat

  if(id_column != "symbol"){
    colnames(mat)[colnames(mat)==id_column] <- "symbol"
  }

  if (sum(duplicated(mat$symbol)) == 0){
    print("No duplicate symbols")
    return(input)
  }

  print(paste("Starting with gene expression matrix containing", nrow(mat), "rows."))

  #Make frequency table
  id_table <- as.data.frame(table(mat$symbol))

  #Identify duplicate genes
  dups <- id_table$Var1[id_table$Freq > 1]
  stopifnot(length(dups) == length(unique(dups)))
  print(paste("Number of genes with duplicate names:", length(dups)))

  #Set aside rows with unique gene names
  nodup_df <- mat[!mat$symbol %in% dups,]

  #Set aside rows with duplicate ids
  dup_df <- mat[mat$symbol %in% dups,]
  stopifnot(nrow(nodup_df) + nrow(dup_df) == nrow(mat))

  #Sort by recurring id
  dup_df <- dup_df[order(dup_df$symbol),]

  print(paste("Number of rows with duplicate gene ids:", nrow(dup_df)))

  #Mean expression fpkm of genes with the same symbol
  mean_exps <- matrix(ncol = ncol(dup_df)-1, nrow=0) #Empty matrix, -1 gene symbol column

  for (i in 1:length(unique(dup_df$symbol))){
    #Subset rows with same symbol, discard symbol column, then apply aggregate function
    exp <- f(as.matrix(dup_df[dup_df$symbol==unique(dup_df$symbol)[i], -1]))
    mean_exps <- rbind(mean_exps, exp)
  }
  stopifnot(nrow(mean_exps) == length(unique(dup_df$symbol)))

  rownames(mean_exps) <- unique(dup_df$symbol)
  mean_exps <- as.data.frame(mean_exps) %>% rownames_to_column("symbol")

  dedupped_df <- rbind(mean_exps, nodup_df)
  dedupped_df <- dedupped_df[order(dedupped_df$symbol),]

  stopifnot(length(unique(dedupped_df$symbol))==length(dedupped_df$symbol)) #All symbols should not be unique
  stopifnot(nrow(mat) - #starting number
              nrow(dup_df) + #rows with duplicate genes...
              length(dups) == #...which condense down into this many unique genes...
              nrow(dedupped_df)) #...should equal the number of rows in the final matrix

  print(paste("Number of genes after applying", substitute(f),  "to duplicate ids:", nrow(dedupped_df)))

  #For estimate, the column with identifiers HAS to be called GeneSymbol or EntrezGeneID
  colnames(dedupped_df)[colnames(dedupped_df)=="symbol"] <- final_gene_symbol_colname

  return(dedupped_df)
}

geneEx = df %>% select(-ensembl_gene_id)
geneEx = summarize_expression_duplicate_ids(mat = geneEx, id_column = "gene_name")
rownames(geneEx) = NULL
geneEx = as.matrix(column_to_rownames(geneEx, "GeneSymbol"))

#### Comparison setup ----

#In order of interest, i.e., which ones we want to have first
comparisons = c("inv_vs_nonprbc", "inv_vs_rest", "prbc_vs_nonprbc", #Most interesting
                'prbc_vs_rest', "lac_vs_rest", "nonprbc_vs_rest", #Remaining one vs rest comparisons
                "inv_vs_prbc", "inv_vs_lac","lac_vs_nonprbc", "lac_vs_prbc") #Remaining pairwise comparisons

sampledata <- as.data.frame(colData(dds))

print("") #Empty line

#### Normalization and voom transformation ----

for (i in 1:length(comparisons)){

  comp = comparisons[i]
  print(paste("Current comparison:", comp))

  #Assemble dge
  y = sampledata[,c("batch", "PAM50", "tumor_purity", comp)]
  dge <- DGEList(counts = geneEx, samples = y, group = y[,comp])

  print(paste("Number of samples", ncol(dge)))
  print(c("Levels in comparison group:", paste(levels(dge$samples$group), collapse=", ")))

  #For pairwise comparisons, we have two options: perform normalization before or after removing samples that aren't part of the comparison
  #To increase statistical sensitivity, we remove them before TMM/voom
  if (str_detect(comp, "vs_rest") == F){
    #y = y[y[,comp] != "other",]
    print(paste("Removing other study groups from comparison", comp))
    dge = dge[,colnames(dge)[dge$samples$group != "other"]]
    print(paste("Number of samples after removing other", ncol(dge)))
    print(c("Levels in comparison group after removing other:", paste(levels(dge$samples$group), collapse=", ")))

  }

  # Store which samples are used for later
  comparison_samples = dge$samples %>% rownames_to_column("sample_name") %>%
    left_join(., select(sampledata, sample_name, patient_ref, study_group, sample_ref, sample_id),
              by = "sample_name")

  print("Performing TMM normalization")
  dge <- calcNormFactors(dge, method = "TMM")
  #print(head(dge$samples))

  print("Performing voom transformation")
  m <- model.matrix(~batch + PAM50 + tumor_purity + group, dge$samples)
  v <- voom(dge, design=m, plot=T, save.plot=T)

  # Save objects
  dir.create(here("results","flexgsea","limma_input"), showWarnings = F, recursive = T)
  saveRDS(list(model = m, voom = v),
          file=here("results", "flexgsea", "limma_input", paste0(comp,"_input_flimma.Rds")))

  #Results dir
  dir.create(here("results","flexgsea","limma_results"), showWarnings = F, recursive = T)

  #### Flexgsea ----

  for (i in 1:length(gene_signatures)){

    sig_name = names(gene_signatures)[i]
    signature = gene_signatures[[i]]

    print(paste("Starting FlexGSEA for comparison", comp, "with gene signature", sig_name))

    start = Sys.time()

    flexres = flexgsea(x = v, y = m,
                       gene.sets = signature,
                       gene.score.fn = flexgsea_limma,
                       es.fn=flexgsea_weighted_ks,
                       sig.fun=flexgsea_calc_sig,
                       nperm=nperm, parallel = T)
    end = Sys.time()
    print(end - start)

    parameters = enframe(
      c(
        signature_name = sig_name, signature_files = signature_files$file[i],
        comparison_name = comp,
        time_completion = format(end, "%a %b %X %Y"), time_start = format(end, "%a %b %X %Y"),
        time_diff_min = (end-start)/60,
        number_permutations = nperm
      ),
      "parameter", "value")


    saveRDS(list(parameters = parameters, sampledata = comparison_samples, flexgsea_results = flexres),
            file=here("results", "flexgsea", "limma_results",
                      paste(comp, sig_name, "results_flimma.Rds", sep="_")))

    print("") #Empty line

    print(paste("Showing results for coefficient:", names(flexres$table)[length(flexres$table)]))
    sigflexres = flexres$table[[length(flexres$table)]] %>% arrange(fdr) %>% filter(fdr < 0.25)

    print(paste(nrow(sigflexres), "gene sets with fdr < 0.25",
                "in comparison", comp, "with gene signature", sig_name))

    if (nrow(sigflexres) > 0){
      print(sigflexres)
    }

    print("") #Empty row
  }
}

sink()

