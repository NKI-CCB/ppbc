#Calc_x functions by Evert Bosdriesz
#plot_enrichment & fisher_pathways wrapper by KSM

#' Calculate enrichment.
#'
#' Calculate enrichment of 'pathway' genes in 'sample' compared to 'population'
#' using hypergeometric test.
#'
#' @param sample - character vector, IDs the selection of genes (from the
#'   population that you want to test for representation in your pathway)
#' @param pathway - character vector, IDs pathway genes
#' @param population - character vector, many genes that represent the
#'   population
#' @export

calc_enrichment <- function(sample, pathway, population) {
  Lpath <- length(unique(pathway[!is.na(pathway)]))
  Lpop <- length(unique(population[!is.na(population)]))
  Ls <- length(unique(sample[!is.na(sample)]))

  # number of pathway genes in sample
  Mps <- as.character(intersect(pathway, sample))
  Lps <- length(Mps)
  Lppop <- length(intersect(pathway, population))
#
#   if ( Lppop < Lpath) {
#     warning('Some pathway members are not in population:\n')
#     # cat(setdiff(pathway, population))
#     }

  pval <- phyper(Lps - 1, Lpath, Lpop - Lpath, Ls, lower.tail = FALSE)
  res <- list(Lps, Ls,  Lpath, Lpop, (Lps/Ls)/(Lppop/Lpop), pval)
  if (length(Mps) == 0) {
    res <- c(res, "None")
  } else{
    res <- c(res, str_c(Mps, collapse = ";"))
  }


  names(res) <- c(
    'path_in_sample', 'tot_in_sample',
    'path_in_pop', 'tot_in_pop',
    'fold_enrichment', 'pval',
    'members_in_sample'
    )

  return(res)
}

#' Calculate enrichment of a collection of gene sets
#'
#' Calculate enrichment of genes in "sample" compared to 'population', for the
#' genesets in "genesets" using hypergeometric test.
#'
#' @param sample - character vector, IDs the selection of genes (from the
#'   population that you want to test for representation in your pathway)
#' @param population - character vector, many genes that represent the
#'   population
#' @param genesets - List. each element is a geneset against which the sample is
#' tested.
#' @param nmin - Minimal number of genes in geneset for the geneset to be used.
#' Default = 10
#' @param nmax - Maximal number of genes in geneset for the geneset to be used.
#' Default = Inf
#' @param fdr_cutoff - Only return results with an fdr < fdr_cutoff
#' @return data_frame object.
#' @export
calc_enrichment_sets <- function(sample, population, genesets,
                                 nmin = 10,
                                 nmax = Inf,
                                 fdr_cutoff = Inf){
  # Discard genesets with < nmin or > nmax samples to aviod multiple testing problems
  ids_use <- names(genesets[sapply(genesets, length) >= nmin & sapply(genesets, length) <= nmax])

  result <- genesets[ids_use] %>%
    purrr::map(function(gs) calc_enrichment(
      sample = sample, pathway = gs, population = population)
      ) %>%
    bind_rows(.id = "pathway") %>%
    dplyr::mutate(fdr = p.adjust(.$pval, method = 'BH')) %>%
    dplyr::filter(fdr < fdr_cutoff) %>%
    dplyr::arrange(pval) %>%
    dplyr::select(-members_in_sample, members_in_sample)

    return(result)

}

#Wrapper function for performing the above on a list of gene signatures
#' @param sig_genes - A list of gene ids that match the ids in the list of pathways
#' @param background_genes - All genes in the dataset
#' @param list_signatures - A list of gene sets that can be provided to calc_enrichment_sets
#' @param verbose - Print number of hits per signature
#' @param collapse_rows - Whether to return a single data frame instead of a list
#' @param fdr_cutoff - Only return results with an fdr < fdr_cutoff
#' @return data_frame object.
#' @export
fisher_pathways <- function(sig_genes, background_genes, list_signatures, fdr_cutoff= Inf, verbose=T, collapse_rows=F){
  res = list()
  for(i in 1:length(list_signatures)){
    this_sig = names(list_signatures)[i]
    #print(paste("Fisher's exact test for", this_sig))
    this_res = calc_enrichment_sets(
      sample = sig_genes,
      genesets  = list_signatures[[i]],
      population = background_genes,
      fdr_cutoff = fdr_cutoff)
    this_res = this_res %>% mutate(collection = this_sig)
    this_res = this_res %>% select(pathway, fdr, collection, everything())
    if(verbose){
      print(paste(nrow(this_res[this_res$fdr < 0.05,]), "significant (FDR < 0.05) pathways in", this_sig))
    }
    res[[this_sig]] = this_res
  }
  if(collapse_rows){
    res = res %>% bind_rows(.id="pathway")
  }
  return(res)
}

#' Dotplot of results for calc_enrichment_sets
#'
#' @param enrich_res - data frame produced by calc_enrichment_sets
#' @param fdr - only plot results below this threshold
#' @param max_nchar_path - for extremely long pathway names, truncate after n characters
#' @param max_path - max number of pathways to plot
#' @param title - the plot title
#' @return ggplot object
#' @export
plot_enrichment <- function(enrich_res, fdr = 0.05, max_nchar_path = 40, max_path = 20, title=""){

  enrich_res = enrich_res %>% arrange(fdr)
  enrich_res = head(enrich_res, max_path)

  enrich_res = enrich_res %>%
    filter(fdr < !!fdr) %>%
    mutate(path = str_remove(pathway, "GO_"),
           Count = path_in_sample) %>%
    mutate(path = if_else(nchar(path) > !!max_nchar_path,
                          str_sub(path, start = 1, end = !!max_nchar_path), path)) %>%
    arrange(Count)
  #Fold enrichment = (N sig genes in pathway/N sig genes)/(Length pathway/All genes)
  enrich_res %>%
    ggplot(aes(x=Count, y = factor(enrich_res$path, levels=enrich_res$path), color = fdr, size = fold_enrichment)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE)) +
    ylab("Path") + ggtitle(title)
}


#' Calculate enrichment of a collection of gene sets for multiple clusters
#'
#' Calculate enrichment of genes in clusters compared to 'population', for the
#' genesets in "genesets" using hypergeometric test.
#'
#' @param sample_clustering - dataframe with 2 columns. First column is feature
#' id, second column is cluster membership of that feature
#' @param population - character vector, many genes that represent the
#'   population
#' @param genesets - List. each element is a geneset against which the sample is
#' tested.
#' @param nmin - Minimal number of genes in geneset for the geneset to be used.
#' Default = 10
#' @param nmax - Maximal number of genes in geneset for the geneset to be used.
#' Default = Inf
#' @param fdr_cutoff - Only return results with an fdr < fdr_cutoff
#' @return data_frame object.
#' @export
calc_enrichment_clustering <- function(sample_clustering,
                                       population,
                                       genesets,
                                       nmin = 10,
                                       nmax = Inf,
                                       fdr_cutoff = Inf){

  stopifnot(ncol(sample_clustering) == 2)

  colnames(sample_clustering) <- c('id', 'cluster')
  clusters <- unique(sample_clustering$cluster)

  # Check if all feature ids in population
  stopifnot(length(setdiff(sample_clustering$id, population)) == 0)

  df_lst <- vector("list", length(clusters))
  names(df_lst) <- clusters

  for (c in clusters) {
    sample <- sample_clustering %>%
      dplyr::filter(cluster == c) %>%
      dplyr::pull("id") %>%
      unique()

    df_lst[[c]] <- calc_enrichment_sets(sample, population, genesets,
                                        nmin = nmin, nmax = nmax,
                                        fdr_cutoff = fdr_cutoff)
    df_lst[[c]]$cluster <- c
  }
  df <- dplyr::bind_rows(df_lst) %>%
    dplyr::mutate(fdr = p.adjust(.$pval, method = 'BH')) %>%
    dplyr::arrange(pval) %>%
    dplyr::select(-members_in_sample, members_in_sample)

  return(df)
}

