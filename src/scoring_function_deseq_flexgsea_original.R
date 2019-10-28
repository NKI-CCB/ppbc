score_genes_deseq <- function(x, y, abs) {
  dds <- DESeqDataSetFromMatrix(t(x), y, formula(paste0('~ ', args$var)))
  dds <- DESeq(dds)
  p <- resultsNames(dds) %>%
    discard(~ . == 'Intercept') %>%
    set_names(., .) %>%
    map_dfc(~ lfcShrink(dds, coef=., type='apeglm',
                        apeMethod='nbinomC')$pvalue) %>%
    data.matrix()
  p[is.na(p)] <- 1.0
  if (abs) {
    1-p
  } else {
    1-p  # FIXME
  }
}