library(DESeq2)
library(tidyverse)
#library(doParallel) #NOPE
library(doMC)
library(flexgsea) #https://github.com/NKI-CCB/flexgsea-r
#library(devtools)
#install_github("NKI-CCB/flexgsea-r", "fix_running_es_at")
library(here)

#BEFORE RUNNING, set environment variable to prevent unwanted BLAS behavior
#export OMP_NUM_THREADS=1 #Run in terminal prior to invoking script, then run script directly from terminal

# Consider: Addiing 'running_es_pos', 'running_es_neg', 'running_es_at' to return_values in flexgsea for producing mountain plots

#Whether to overwrite existing results file
overwrite=F

dir.create(here("results", "flexgsea", "logs"), showWarnings = F, recursive = T)

logfile <- file(here("results", "flexgsea", "logs", "flexgsea_deseq.txt"),)
sink(logfile ,type = "output")
#sink(logfile, type = "message") #Doesn't work properly

#Set number of permutations for flexgsea
nperm = 1000 #5 for testing, 1000 for real thing
print(paste("Number of permutations:", nperm))

#Set output dirs

# Input before running function
DIR = here("results","flexgsea", "deseq")
input_DIR = file.path(DIR,"input")
dir.create(input_DIR, showWarnings = F, recursive = T)

#Results dir
results_DIR = file.path(DIR,"results")
dir.create(results_DIR, showWarnings = F, recursive = T)

#registerDoParallel(cores=32) #doParallel version
doMC::registerDoMC(32)

source(here("src", "deseq_report_functions.R"))
#source(here("src", "scoring_function_deseq_flexgsea.R")) We use a lambda function instead

####Load data----

#Pre-filtered, with all necessary column groups
dds = readRDS(here("data/Rds/08_dds_inv_vs_rest_standard.Rds"))
#testing
#test.dds = head(dds,1000)
gx_annot <- read_tsv(here("data/metadata/01_tx_annot.tsv"))
gx_annot = gx_annot %>% select(ensembl_gene_id = gene_id, gene_name, gene_type, description = gene_description) %>% distinct()

#Load gene signatures

gene_signatures = list(
  go_c5 = read_gmt(file = here("data","external","gmt","c5.all.v7.0.symbols.gmt")),# Gene ontology BP
  hallmark = read_gmt(file = here("data","external","gmt","h.all.v7.0.symbols.gmt")), # 50 gene overview
  canonpath_c2 = read_gmt(file = here("data","external","gmt","c2.cp.v7.0.symbols.gmt")), #Canonical pathways
  cgp_c2 = read_gmt(file = here("data", "external", "gmt", "c2.cgp.v7.0.symbols.gmt")) #Chemical and genetic perturbations
  #oncogene_c6 = read_gmt(file = here("data","external","gmt","c6.all.v6.2.symbols.gmt")), #Oncogenic gene signatures
  #immune_c7 = read_gmt(file = here("data","external","gmt","c7.all.v6.2.symbols.gmt")), #Immunological signatures
  #motif_tf_mir_c3 = read_gmt(file = here("data","external", "gmt","c3.all.v6.2.symbols.gmt")), #Targets of TFs and miRs
  #comp_canc_c4 = read_gmt(file = here("data","external","gmt","c4.all.v6.2.symbols.gmt")) #Computational gene sets from cancer gene neighborhoods and cancer modules
)

signature_files = enframe(c(
  hallmark = here("data","external","gmt","h.all.v6.2.symbols.gmt"), # 50 gene overview
  canonpath_c2 =  here("data","external","gmt","c2.cp.v6.2.symbols.gmt"), #Canonical pathways
  cgp_c2 = here("data", "external", "gmt", "c2.cgp.v7.0.symbols.gmt"), #Chemical and genetic perturbations
  go_c5 = here("data","external","gmt","c5.all.v6.2.symbols.gmt")),
  "signature_name", "file")

names(gene_signatures) %>% head()

#Define whether each gene signature should be run in absolute mode or not
#From Tycho:
#Hallmark performs better when not absolute, and cp better when absolute
#CGP is obviously directional
#GO terms are sometimes directional and sometimes not, go with absolute

signature_abs = tibble(
  signature = names(gene_signatures),
  abs = c(F, T, F, T))

print("Absolute GSEA or not, by signature")
print(signature_abs)


#### Comparison overview ----

print("Pairwise comparisons:")
print(table(dds$study_group))
print("One vs rest comparisons")
print(table(dds$inv_vs_rest, dds$study_group))
print(table(dds$prbc_vs_rest, dds$study_group))
print(table(dds$lac_vs_rest, dds$study_group))
print(table(dds$nonprbc_vs_rest, dds$study_group))


### Convert to gene symbol count matrix ----

#Extract raw counts
counts <- assay(dds)
#Testing
#counts.test = assay(test.dds)

#Gmt files never have Ensembl gene ids, so we'll need to convert to gene symbols.

df <- right_join(gx_annot[,c("ensembl_gene_id","gene_name")],rownames_to_column(as.data.frame(counts),"ensembl_gene_id"),
                 by="ensembl_gene_id")

print(paste("Number of ensembl IDs with no corresponding gene symbol:", sum(is.na(df$gene_name))))
df <- df[!is.na(df$gene_name),]

#Testing
#test.df <- right_join(gx_annot[,c("ensembl_gene_id","gene_name")],rownames_to_column(as.data.frame(counts.test),"ensembl_gene_id"),
#                      by="ensembl_gene_id")
#test.df <- test.df[!is.na(test.df$gene_name),]

#De-duplicate genes
geneEx = df %>% select(-ensembl_gene_id)

#When dealing with duplicated gene names, take the mean expression
geneEx = summarize_expression_duplicate_ids(mat = geneEx, id_column = "gene_name", f=colMeans)
saveRDS(geneEx, file.path(input_DIR, "geneSymbol_countmatrix.Rds"))

#Testing
#test.geneEx <- test.df %>% select(-ensembl_gene_id)
#test.geneEx = summarize_expression_duplicate_ids(mat = geneEx, id_column = "gene_name", f=colMeans) #No dups
#test.geneEx <- column_to_rownames(test.geneEx, "gene_name")

#### Comparison setup ----

#In order of interest, i.e., which ones we want to have first
sampledata <- as.data.frame(colData(dds))

#testing
#test.sampledata <- as.data.frame(colData(test.dds))

#In order of interest, i.e., which ones we want to have first
comparisons = c("inv_vs_rest",  #Most interesting
                'prbc_vs_rest',
                "study_group",
                "lac_vs_rest",
                "nonprbc_vs_rest"
                )


print("") #Empty line

#### Lambda score genes function ----
make_score_genes <- function(coef_of_interest="study_group") {
  score_genes_deseq <- function(x, y, abs) {

    require(DESeq2)
    require(tidyverse)

    stopifnot(coef_of_interest %in% colnames(y))

    #Required for permutation test
    rownames(x) <- rownames(y)

    genes = colnames(x)
    #Create dds with formula derived from column names
    dds <- DESeqDataSetFromMatrix(t(x), y, as.formula(paste("~", paste(colnames(y), collapse="+"))))
    #print(design(dds))
    dds <- DESeq(dds)

    #Extract result names
    coefficients = resultsNames(dds) %>%
      discard(~ . == 'Intercept') %>%
      set_names(., .)

    #Subset down to column(s) of interest
    if (!is.null(coef_of_interest)){

      #Retrieve the levels from each of the columns of interest
      ourgroups = as.character(unlist(apply(y[,coef_of_interest, drop=F], 2, unique)))

      #Subset results names by pattern matching to levels in columns of interest
      #Mapping returns multiple non-unique matches as levels will appear more than once in the design matrix
      #Use ~keep instead of ~detect
      #map_chr expects a single match and will crash when multiple hits are found
      #Solution: Use map instead of map_chr and unlist
      coefficients =
        unique(
          unlist(map(ourgroups, ~ keep(coefficients, str_detect, .x)))
        )
      #coefficients
    } else {
      coefficients = tail(coefficients, 1) #Last one in group
    }

    #Shrink coefficients and extract pval and l2fc
    res <- coefficients %>% map_dfc(~ as.data.frame(lfcShrink(dds, coef=., type='apeglm',
                                                              apeMethod='nbinomC'))[,c("pvalue","log2FoldChange")]) #Include l2fc so we can adjust by sign
    rownames(res) = genes
    #Separate p val
    p = res %>% select(contains("pval"))
    colnames(p) = coefficients
    p = data.matrix(p)
    p[is.na(p)] <- 1.0

    #Create data frame indicating whether l2fc is above or below 0
    l2fc = res %>% select(contains("log2FoldChange"))
    colnames(l2fc) = coefficients
    l2fc = data.matrix(l2fc)
    sign = ifelse(l2fc > 0, 1, -1)

    if (abs) {
      return(1-p)
    } else {
      return((1-p) * sign )} #Multiply pval by -1 if l2fc < 0
  }}


#### Flexgsea ----

x = t(geneEx)
mode(x)="integer" #Crashes if numeric
saveRDS(geneEx, file.path(input_DIR, "x.Rds"))

#testing
#test.x <- t(test.geneEx)
#mode(test.x) <- "integer"
#saveRDS(test.x, file.path(input_DIR, "test.x.Rds"))


for (i in 1:length(comparisons)){
  comp = comparisons[i]

  for (i in 1:length(gene_signatures)){

    sig_name = names(gene_signatures)[i]
    signature = gene_signatures[[i]]

    y = sampledata %>% select(batch, PAM50, !!comp)
    #testing
    #test.y=test.sampledata %>% select(batch, PAM50, !!comp)

    print(paste("Starting FlexGSEA for comparison", comp, "with gene signature", sig_name))

    #Test whether analysis has been run before
    result_file = file.path(results_DIR, paste(comp, sig_name, "results_flexdeseq.Rds", sep="_"))

    if(file.exists(result_file) & overwrite == F){
      print("Results file exists and overwrite is false, skipping to next")
      next
    }

    if(file.exists(result_file) & overwrite == T){
      print("Results file exists and overwrite is true, replacing pre-existing results file")
    }

    start = Sys.time()

    #To save time, we design the scoring function to return pairwise comparisons simultaneously
    #For one vs rest, selects the last column only (which corresponds to one vs rest)
    if (comp == "study_group"){
      this.comp = "study_group"
    } else {
      this.comp = NULL
    }


    #Should this run be absolute?
    this_abs = as.logical(signature_abs[signature_abs$signature == sig_name,"abs"])

    flexres = flexgsea(x = x, y = y,
                       gene.sets = signature,
                       gene.score.fn = make_score_genes(coef_of_interest=this.comp),
                       es.fn=flexgsea_weighted_ks,
                       sig.fun=flexgsea_calc_sig,
                       nperm=nperm, parallel = T,
                       block.size = 1,
                       abs = this_abs,
                       return_values = c("gene_name","leading_edge",'running_es_pos','running_es_neg')
                       )
    #'running_es_at' = Error in res$running_es_at[[i]] <- list() : replacement has length zero
    #test.flexres = flexgsea(x = test.x, y = test.y,
     #                  gene.sets = signature,
      #                 gene.score.fn = make_score_genes(coef_of_interest=this.comp),
       #                es.fn=flexgsea_weighted_ks,
        #               sig.fun=flexgsea_calc_sig,
         #              nperm=nperm, parallel = T,
          #             block.size = 1,
           #            abs = this_abs,
            #           return_values = c("gene_name","leading_edge",'running_es_pos','running_es_neg')
    #)

    end = Sys.time()
    print(end - start)

    parameters = enframe(
      c(
        signature_name = sig_name, signature_files = signature_files$file[i],
        comparison_name = comp,
        time_completion = format(end, "%a %b %X %Y"), time_start = format(end, "%a %b %X %Y"),
        time_diff_min = (end-start)/60,
        number_permutations = nperm,
        abs = this_abs
      ),
      "parameter", "value")


    saveRDS(list(parameters = parameters, flexgsea_results = flexres),
            file=result_file)


    #Testing
    #saveRDS(list(parameters = parameters, flexgsea_results = test.flexres),
    #        file=file.path(results_DIR,
    #                       paste(comp, sig_name, "testresults_flexdeseq.Rds", sep="_")))

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
#sink() #Close both warning and output connections
