
#' Score genes with DESeq
#'
#' Uses  p value to rank the genes (after shrinking the fold changes with apeglm)
#' P val is used instead of padj because adjustment methods force a lot of values to 1, causing loss of information
#' If abs == F, pval is multiplied by -1 if log2fc < 0
#'
#' @param x - An expression matrix formulated for flexgsea
#'   (samples in rows and genes in columns)
#' @param y - A response variable matrix for flexgsea, i.e. the columns in
#'   colData that make up your dds design
#' @param abs - Whether to detect the directionality of the gene
#'   by multiplying by the sign of the l2fc
#' @param coef_of_interest - One or more columns in y. Some design formulas may contain
#'   elements that should be excluded from pathway analysis (i.e. batch). You can select
#'   variables of interest by specifying the column name.
#'   Because we can't pass arguments to flexgsea, default is the last coefficient in the design formula
#' @export


#Why not rank with the Wald stat? Because it's not reported when apeglm is the shrinkage method (and no apparent way to retrieve it)

#Flexgsea can't take parameters with the scoring function
#We therefore have two versions: One that just takes the last coefficient, and one which takes all levels of the comparison of interest


score_genes_deseq <- function(x, y, abs, coef_of_interest) {

  require(DESeq2)
  require(tidyverse)

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
  if (coef_of_interest != "last"){

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
  }




#Testing

#dds <- readRDS(file = here("data/Rds/05_dds_PAM50_batch.Rds")) #All comparisons with study_group are relevant
#tiny.dds = head(dds)
#design(tiny.dds)

#x = t(assay(tiny.dds))
#y <- as.data.frame(colData(tiny.dds)) %>% select(batch, PAM50, study_group) #Pairwise comparisons
#resultsNames(DESeq(tiny.dds))
#test = score_genes_deseq(x,y, abs=F)
#test = score_genes_deseq(x,y, abs=F, coef_of_interest = "study_group")
#test
#test = score_genes_deseq(x,y, abs=F, coef_of_interest = c("PAM50","study_group"))
#test

#dds = readRDS(here("data/Rds/08_dds_inv_vs_rest_standard.Rds")) #Has only one group of interest (inv vs rest)
#tiny.dds = head(dds)
#design(tiny.dds)

#x = t(assay(tiny.dds))
#y <- as.data.frame(colData(tiny.dds)) %>% select(batch, PAM50, inv_vs_rest) #One vs rest
#test = score_genes_deseq(x,y, abs=F, coef_of_interest = "inv_vs_rest")
#test
