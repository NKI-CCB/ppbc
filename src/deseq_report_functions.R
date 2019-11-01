#### Enrichment analysis functions
#Hypergeom tests by Evert Bosdriesz
source(here("src", "enrichment-analysis-functions.R"))


#### Standard report functions ####

annotate_results <- function(results_object, anno_df = gx_annot, immune_list = immune_gene_list,
                             mark_immune = T, mark_outliers = T){

  #Takes a res object from deseq, merges it with a gene annotation file, and an optional immune list to create a data frame sorted by adjusted p value
  #Type is an aggregate of immunoglobulin genes from gene type and genes with a known immune function from ImmPort
  require(tidyverse)

  anno_res <- results_object %>% as.data.frame() %>% rownames_to_column("ensembl_gene_id") %>%
    right_join(anno_df,., by = "ensembl_gene_id") %>% arrange(padj)

  if (mark_immune==T){ #Shows genes that defined as immune based on external list
    anno_res = anno_res %>% mutate(ImmPort_gene = if_else(
      ensembl_gene_id %in% immune_list$ensembl, T, F
    ))

    anno_res = anno_res %>%
      mutate(Type=case_when( #Order is hierarchical
        str_detect(string=gene_type,pattern="^IG") ~ "immunoglobulin",
        str_detect(string=gene_type,pattern="^TR") ~ "T cell receptor",
        ImmPort_gene == T ~ "immune protein coding",
        gene_type == "protein_coding" ~ "protein coding",
        TRUE ~ "other noncoding"))

    anno_res = anno_res %>% select(gene_name,gene_type,Type, everything())
  }

  #DESeq sets both the pval and padj to NA for genes above a Cooks distance threshold
  #Cook’s distance is a measure of how much a single sample is influencing the fitted coefficients for a gene,
  #and a large value of Cook’s distance is intended to indicate an outlier count.
  #If a row is filtered by automatic independent filtering, for having a low mean normalized count, then only the adjusted p value will be set to NA.
  #Keep track of both possibilities

  if (mark_outliers == T){
    anno_res = anno_res %>%
      mutate(cooks_outlier = if_else(baseMean > 0 & is.na(pvalue), T, F)) %>%
      mutate(indepfilt = if_else(!is.na(pvalue) & is.na(padj), T, F))
  }


  return(anno_res)
}


significance <- function(annotated_results, pthresh = 0.05, absl2fc = 0.5){
  require(tidyverse)
  sig = annotated_results %>% filter(padj < !!pthresh & (abs(log2FoldChange) > !!absl2fc))
  return(sig)
}

cooks_cutoff <- function(dds){
  #This is how DESeq2 decides the Cooks cutoff by default
  m <- nrow(attr(dds, "dispModelMatrix"))
  p <- ncol(attr(dds, "dispModelMatrix"))
  defaultCutoff <- qf(0.99, p, m - p)
  return(defaultCutoff)
}

vp_cooks_outliers <- function(annotated_results, dds, cook_threshold=25, fc_threshold=0.5, title=""){
  #Visualize a volcano plot of genes which have been set to NA due to Cook's distance
  left_join(annotated_results,enframe(mcols(dds)$maxCooks, "ensembl_gene_id", "max_cooks"), by = "ensembl_gene_id") %>%
    mutate(label = if_else(max_cooks > !!cook_threshold | (cooks_outlier == T & abs(log2FoldChange) > !!fc_threshold),
                           gene_name, NULL)) %>%
    ggplot(aes(x=log2FoldChange, y = max_cooks, color = -log10(padj),
               shape=cooks_outlier, size=log10(baseMean), label=label)) +
    geom_point() + ggrepel::geom_label_repel(show.legend = F) +
    ggtitle(paste("Cooks rejects:", title)) +
    geom_hline(yintercept = cook_threshold, color = "black", linetype="dashed")
}

cooks_replacement <- function(dds){
  #When there are 7 or more replicates for a given sample, the DESeq function will automatically replace counts
  #that have Cooks distance greater than the default threshold
  #with the trimmed mean over all samples, scaled up by the size factor or normalization factor for that sample.
  #The original counts are preserved in counts(dds)
  cooks_replaced = mcols(dds)$replace %>% enframe("ensembl_gene_id", "outlier_replaced") %>% filter(outlier_replaced==T) %>%
    left_join(.,select(gx_annot, "ensembl_gene_id", "gene_name", "gene_type"), by = "ensembl_gene_id")
  return(cooks_replaced)
}

summarize_expression_duplicate_ids <- function(mat, id_column=NULL, f=colMeans, verbose=F){
  require(dplyr)

  #Easiest way to write functions with dplyr is to standardize the column name

  input = mat #Save input before doing anything to it in case there are no duplicates

  if(is.null(id_column)==T){
    mat = as.data.frame(mat)
    mat = rownames_to_column(mat, "symbol")
    id_column = "symbol"
  }

  if(id_column != "symbol"){
    colnames(mat)[colnames(mat)==id_column] <- "symbol"
  }

  if (sum(duplicated(mat$symbol)) == 0){
    print("No duplicate symbols")
    return(input)
  }

  if(verbose==T){
    print(paste("Starting with gene expression matrix containing", nrow(mat), "rows."))
    }

  #Make frequency table
  id_table <- as.data.frame(table(mat$symbol))

  #Identify duplicate genes
  dups <- id_table$Var1[id_table$Freq > 1]
  stopifnot(length(dups) == length(unique(dups)))
  if(verbose == T){
    print(paste("Number of genes with duplicate names:", length(dups)))
    }

  #Set aside rows with unique gene names
  nodup_df <- mat[!mat$symbol %in% dups,]

  #Set aside rows with duplicate ids
  dup_df <- mat[mat$symbol %in% dups,]
  stopifnot(nrow(nodup_df) + nrow(dup_df) == nrow(mat))

  #Sort by recurring id
  dup_df <- dup_df[order(dup_df$symbol),]


  if(verbose==T){
    print(paste("Number of rows with duplicate gene ids:", nrow(dup_df)))
  }

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

  if (verbose == T){
    print(paste("Number of genes after applying", substitute(f),  "to duplicate ids:", nrow(dedupped_df)))
  }


  rownames(dedupped_df) = NULL #Required for column_to_rownames
  dedupped_df = column_to_rownames(dedupped_df, "symbol")

  return(dedupped_df)
}

volcano_plot <- function(annotated_results, title="", path_save_fig = NULL, pthresh = 0.05, absl2fc=0.5,
                         shape_col=NULL, dedup_ids = T, colors = gene_colors){

  require(tidyverse)
  require(ggrepel)
  require(here)

  volcdf <- annotated_results %>% mutate(NegLog10FDR = -log10(padj),
                                         Significant = if_else(
                                           (padj < !!pthresh & abs(log2FoldChange) > !!absl2fc),
                                           TRUE, FALSE),
                                         #Type=case_when( #Redundant
                                          # gene_type == "protein_coding" ~ "protein coding",
                                           #str_detect(string=gene_type,pattern="IG") ~ "IG gene",
                                           #TRUE ~ "other noncoding"),
                                         Color = if_else(
                                           Significant == T, Type, "n.s."
                                         )) %>%
    filter(!is.na(Significant))

  if (dedup_ids == T){
    print("Summarizing -log10FDR and log2FoldChange by median for duplicate gene names in volcano plot")
    volcdf = volcdf %>% group_by(gene_name, Significant, Type, Color) %>%
      summarise(log2FoldChange = median(log2FoldChange), NegLog10FDR = median(NegLog10FDR)) %>%
      arrange(desc(NegLog10FDR))
  }

  colors = c(colors,"n.s."="gray")

  volcplot <-ggplot(volcdf, aes(x=log2FoldChange,y=NegLog10FDR, color=Color)) +
    geom_point() +
    scale_color_manual(values=colors) +
    geom_hline(yintercept = -log10(pthresh), colour="darkred", linetype="dashed") +
    geom_vline(xintercept = absl2fc, color="darkred", linetype="dashed") +
    geom_vline(xintercept = -absl2fc, color="darkred", linetype="dashed") +
    theme_bw()

  volcplot <-
    volcplot +
    geom_text_repel(data=head(volcdf[volcdf$Significant==T,], 10), aes(label=gene_name), show.legend = F) +
    ggtitle(title) +
    labs(y="-log10(FDR)")

  if (is.null(path_save_fig)==F){
    suppressMessages(ggsave(filename=path_save_fig, plot = volcplot))
  }

  if (is.null(shape_col)==F){
    volcplot = volcplot + aes(shape=get(shape_col)) + labs(shape=shape_col)
  }

  return(volcplot)

}

color_grid = function (colours, labels = T, names = T, borders = NULL, cex_label = 1)
{

  #Adapted from scales' show_colors
  #Visualizes colors chosen for each study group

  n <- length(colours)
  ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n/ncol)

  #Also add option to plot group labels
  colornames = c(names(colours), rep(NA, nrow * ncol - length(names(colours))))
  colornames <- matrix(colornames, ncol = ncol, byrow = TRUE)

  #Original scales code
  colours <- c(colours, rep(NA, nrow * ncol - length(colours)))
  colours <- matrix(colours, ncol = ncol, byrow = TRUE)


  old <- par(pty = "s", mar = c(0, 0, 0, 0))
  on.exit(par(old))
  size <- max(dim(colours))
  plot(c(0, size), c(0, -size), type = "n", xlab = "", ylab = "",
       axes = FALSE)
  rect(col(colours) - 1, -row(colours) + 1, col(colours), -row(colours),
       col = colours, border = borders)



  if (labels) {
    text(col(colours) - 0.5, -row(colours) + 0.3, colours,
         cex = cex_label)
  }

  #Plot the group labels
  if (names) {
    text(col(colornames) - 0.5, -row(colornames) + 0.5, colornames,
         cex = cex_label)
  }
}





deseq_heatmap = function(mat, sampledata, sig_results,
                         groups_to_plot=levels(sampledata$study_group),
                         title=NULL,
                         top_colors = study_colors,
                         bottom_colors = pam_colors,
                         row_colors = gene_colors,
                         row_scale = T,
                         maxn_rownames = 50,
                         row_size = 8, col_size = 8,
                         show_col_names = F, ...){
  require(scrime)
  require(ComplexHeatmap)
  require(tidyverse)

  #stopifnot(groups_to_plot %in% levels(sampledata$study_group)) #Now accepts other columns
  stopifnot(identical(colnames(mat), sampledata$sample_name)) #Ensure that a column containing all sample names exists


  if(is.null(title)==T){
    title = if_else(identical(groups_to_plot,levels(sampledata$study_group)),
                    "All groups",
                    paste(groups_to_plot, collapse="vs"))
  }


  #Reduce genes to significant only
  genestoplot = sig_results$gene_name
  mat = mat[rownames(mat) %in% genestoplot, ]

  #Reduce genes and sample data to compared groups only
  sampledata = sampledata %>% filter(study_group %in% groups_to_plot)
  mat = mat[,colnames(mat) %in% sampledata$sample_name]
  stopifnot(identical(colnames(mat), sampledata$sample_name))

  #Toggle row labels depending on size of heatmap
  if (nrow(mat) <= maxn_rownames){
    show_row_names = T
  } else {
    show_row_names = F
  }

  #Row scale settings
  if (row_scale==T){
    mat = scrime::rowScales(mat)
  }
  #Change legend according to whether input is scaled
  if (row_scale==T){
    hlp = list(title="rowscaled vst")
  } else {
    hlp = list(title="vst counts")
  }

  #Ensure input is matrix
  mat = as.matrix(mat)

  #Heatmap annotation
  ann_top = sampledata[,"study_group", drop=F]

  #Top column annotation
  colTop <- HeatmapAnnotation(df=ann_top, which="col",
                              col = list(study_group=top_colors))

  #Bottom column annotation
  ann_bottom = sampledata[,"PAM50", drop=F]
  colBottom <- HeatmapAnnotation(df=ann_bottom, which="col", col = list(PAM50 = bottom_colors))

  # Row annotation
  anno_rows = sig_results %>%
    select(gene_name, Type) %>%
    distinct() %>%
    filter(!duplicated(gene_name)) %>%
    column_to_rownames("gene_name")

  #Essential that the order be the same!
  anno_rows = anno_rows[match(rownames(mat),rownames(anno_rows)), ,drop=F]

  rowAnno = HeatmapAnnotation(df=anno_rows, which="row", col=list(Type = row_colors))

  Heatmap(mat,
          top_annotation = colTop,
          bottom_annotation = colBottom,
          left_annotation = rowAnno,
          heatmap_legend_param = hlp,
          show_row_names = show_row_names,
          show_column_names = show_col_names,
          cluster_rows = T,
          row_names_gp = gpar(fontsize = row_size),
          column_names_gp = gpar(fontsize = col_size),
          column_title = title,
          ...)
}


plot_gene_beehive = function(dds, result_df, groups_to_plot = levels(colData(dds)[,intgroup]),
                             intgroup="study_group", colorby="PAM50", path_save_fig = NULL){

  #Takes a results data frame subsetted down to just the gene of interest
  #For example : plot_gene_beehive(dds, result_df = res_list$LRT[res_list$LRT$ensembl_gene_id == "ENSG00000156738",])

  require(DESeq2)
  require(tidyverse)

  if(nrow(result_df) > 1){
    warning("More than one result provided, plotting first row")
    result_df = head(result_df, 1)
  }

  if(nrow(result_df) < 1){
    stop("Gene ID not in results")
  }

  ensembl_gene_id = result_df$ensembl_gene_id
  gene_name = result_df$gene_name

  stopifnot(groups_to_plot %in% levels(colData(dds)[,intgroup]))

  comparison = if_else(groups_to_plot == levels(colData(dds)[,intgroup]),
                       "all groups", paste(groups_to_plot, collapse = " vs "))


  samples_to_keep = colnames(dds)[as.data.frame(colData(dds)[intgroup])[,1] %in% groups_to_plot]

  dds.trim = dds[,samples_to_keep]

  #Normalizes by size factors
  d <- plotCounts(dds.trim, gene=ensembl_gene_id, intgroup=c(intgroup,colorby), normalize=T, returnData=TRUE)

  title = paste0(paste(gene_name,ensembl_gene_id, sep=":"), ", ", comparison, ", padj: ",
                 formatC(result_df$padj, format = "e", digits = 2),
                 ", log2fc: ", round(result_df$log2FoldChange,2))

  dplot = ggplot(d, aes(x=get(intgroup), y=count)) +
    geom_boxplot(alpha=0, show.legend = F, outlier.shape = NA) +
    geom_point(position=position_jitter(w=0.1,h=0), aes(color=get(colorby))) +
    scale_y_log10() + xlab(intgroup) + labs(color=colorby) + ylab("size-factor normalized counts") +
    ggtitle(title)

  if (is.null(path_save_fig)!=T){
    suppressMessages(ggsave(filename=path_save_fig, plot = volcplot))
  }

  return(dplot)
}

plot_many_beehives = function(dds, result_df, ...){
  plotlist = list()
  for (i in 1 : nrow(result_df)){
    plotlist[[length(plotlist)+1]]= plot_gene_beehive(dds, result_df = result_df[i,], ...)
  }
  return(plotlist)
}


deseq_report = function(results_object, dds, anno_df = gx_annot, mark_immune=T, immune_list=immune_gene_list,
                        pthresh = 0.05, absl2fc = 0.5, variance_stabilized_dds = vsd,
                        title=NULL, list_gene_sets = gene_sets, plot_pathway = T,
                        intgroup="study_group", colorby="PAM50",
                        groups=levels(colData(dds)[,intgroup]),
                        beehive_groups = "comparison",
                        n_beehive_plots=5, verbose=T, ...){
  require(DESeq2)
  require(tidyverse)

  stopifnot(beehive_groups %in% c("all", "comparison"))

  mega = list()

  if(is.null(title==T)){
    title = if_else(identical(sort(groups),
                              sort(levels(colData(dds)[,"study_group"]))), "All groups - Interpret with caution",
                    paste(groups, collapse=" vs "))
  }

  #Merge deseq results object with gene name and type annotation
  if(verbose==T){
    print("Annotating results...")
  }

  df = annotate_results(results_object, mark_immune = mark_immune, immune_list = immune_list)

  mega$annotated_results = df

  #Set significance threshold and save significant hits
  print(paste("Significance threshold is padj:", pthresh, "and abs(log2foldchange):", absl2fc))

  mega$sig_threshold = list(padj = pthresh, absl2fc = absl2fc)

  sig = significance(df, pthresh = pthresh, absl2fc = absl2fc)

  print(paste("Number of unique significant genes:", length(unique(sig$gene_name))))

  mega$significant_genes = sig

  #Draw volcano plot
  if(verbose==T){
    print("Creating volcano plot...")
  }

  vp = volcano_plot(annotated_results = df, title=title, pthresh = pthresh, absl2fc = absl2fc)

  mega$volcano_plot = vp


  #Prepare deduplication
  if(verbose==T){
    print("Converting IDs to symbols and deduplicating repeated symbols...")
  }

  geneEx = rownames_to_column(as.data.frame(assay(variance_stabilized_dds)), "ensembl_gene_id")
  geneEx = right_join(select(anno_df, gene_name, ensembl_gene_id), geneEx, by = "ensembl_gene_id") %>%
    select(-ensembl_gene_id)

  dedup_geneEx = summarize_expression_duplicate_ids(geneEx, id_column = "gene_name")

  #Draw heatmap
  if(verbose==T){
    print("Creating heatmap...")
  }

  hm = deseq_heatmap(mat = dedup_geneEx, sampledata = as.data.frame(colData(variance_stabilized_dds)),
                     sig_results = sig, title = title, groups_to_plot = groups)


  mega$heatmap = hm

  #Draw a selection of beehive plots from top hits
  if(verbose==T){
    print(paste("Creating beehive plots from top", n_beehive_plots, "genes..."))
  }

  if(beehive_groups == "all"){
    beehive_groups = levels(dds$study_group) #Plot all study groups in beehive
  } else {
    beehive_groups = groups #Plot only the ones involved in calculating significance
  }

  beelist = plot_many_beehives(dds = dds, result_df = head(sig, n_beehive_plots), groups_to_plot = beehive_groups)

  mega$beehive_plots = beelist

  #Fisher's exact test for gene signatures
  #sig_genes, background_genes, list_signatures, fdr_cutoff= Inf, verbose=T, collapse_rows=F
  if(!is.null(list_gene_sets)){
    if(verbose==T){print("Fisher's exact test for pathway analysis...")}
    pathways = fisher_pathways(
      sig_genes = sig$gene_name,
      background_genes = df$gene_name,
      list_signatures = gene_sets,
      fdr_cutoff = Inf,
      verbose = verbose,
      collapse_rows=F)
    mega$pathways = pathways
  }

  #Dotplot of pathways
  if(plot_pathway==T & !is.null(list_gene_sets)){
    print("Pathway plots")
    pathway_plots = list()
    for(i in 1:length(pathways)){
      thispath = pathways[[i]]
      if(nrow(thispath[thispath$fdr < 0.05, ]) == 0 ){
        next
      }
      path_plot = plot_enrichment(
        enrich_res = thispath,
        fdr = 0.05,
        max_nchar_path = 40,
        max_path = 20,
        title = paste(title, names(pathways)[[i]])
      )
      pathway_plots[[names(pathways)[i]]] <- path_plot
    }
    mega$pathway_plots <- pathway_plots
  }

  #Outlier reporting
  #Keep track of genes with padj set to NA due to Cooks and independent filtering
  cooks_threshold = cooks_cutoff(dds)
  mega$cooks_threshold = cooks_threshold

  #print(paste("Cook's distance threshold", round(cooks_threshold, 2)))
  print(paste("Genes with padj set to 0 due to exceeding Cooks threshold:", nrow(df[df$cooks_outlier==T,])))
  print(paste("Genes discarded by automatic independent filtering for having a low mean normalized count:",
              nrow(df[df$indepfilt==T,])))

  #Volcano plot of extreme outliers
  mega$volc_outliers <- vp_cooks_outliers(df, dds,
                                          cook_threshold = 25, #Show names for all above default threshold usually overwhelms the plot
                                          fc_threshold = absl2fc,
                                          title=title
                                          )
  #Beehive plots of those which both exceend Cook's threshold AND the l2fc threshold (first 10)
  notable_outliers = df %>% filter(cooks_outlier==T & log2FoldChange > !!absl2fc) %>% pull(ensembl_gene_id)

  outlier_df = df[df$ensembl_gene_id %in% notable_outliers,]
  outlier_df = outlier_df %>% arrange(desc(abs(log2FoldChange)))

  outlier_bees = plot_many_beehives(dds = dds,
                                    result_df = head(outlier_df, 10),
                                    groups_to_plot = beehive_groups)
  mega$outlier_bees = outlier_bees

  return(mega)

}


shrinkRes <- function(dds, contrast, type="apeglm"){
  require(DESeq2)
  require(tidyverse)

  cdf = enframe(colnames(coef(dds)), "coef", "name")
  contrast_name = paste(contrast[1], contrast[2], "vs",contrast[3], sep="_")

  if(!contrast_name %in% cdf$name){
    stop("Contrast name not found within coefficients. Do you need to relevel your design formula?")
  }

  print("Detecting coefficient corresponding to contrast:")

  c = cdf[cdf$name==contrast_name,]
  print(paste("Contrast name:", c$name))
  print(paste("Coefficient:", c$coef))

  print("Retrieving results from deseq dataset with contrast...")
  res = results(dds, contrast=contrast)

  print(paste("Shrinking log2foldchanges via method", type))
  ape <- lfcShrink(dds, coef = c$coef, res = res, type=type)

  return(ape)
}

gene_results = function(list_reports, ensembl_id, pthresh = 0.05, abslogfcthresh = 0.5){

  require(tidyverse)

  #print(paste("Adjusted p value threshold:", pthresh, ", abs(log2FoldChange) threshold:", abslogfcthresh))

  comp = lapply(list_reports,
                function(x) x$annotated_results[x$annotated_results$ensembl_gene_id==ensembl_id,
                                                c("ensembl_gene_id","gene_name","log2FoldChange","padj")]) %>%
    unlist() %>% enframe() %>% separate(name, into = c("comparison", "field"), sep="\\.") %>%
    spread(key=field, value=value) %>%
    mutate(log2FoldChange = as.numeric(log2FoldChange),
           padj = as.numeric(padj),
           padj_cutoff = !!pthresh,
           absLog2FoldChange_cutoff = !!abslogfcthresh)

  comp = comp %>% mutate(sig = (padj < !!pthresh) & (abs(log2FoldChange) > !!abslogfcthresh)) %>%
    select(comparison, sig, padj, log2FoldChange, everything())

  return(comp)
}



blind_vst = function(dds){
  design(dds) = formula(~ 1)
  vsd = vst(dds, blind=T)
  #mat = assay(vsd)
  return(vsd)
  }

gene_report = function(list_reports, ensembl_id, pthresh = 0.05, abslogfcthresh = 0.5){

  require(tidyverse)

  #print(paste("Adjusted p value threshold:", pthresh, ", abs(log2FoldChange) threshold:", abslogfcthresh))

  comp = lapply(list_reports,
                function(x) x[x$ensembl_gene_id==ensembl_id,
                              c("ensembl_gene_id","gene_name","log2FoldChange","padj")]) %>%
    unlist() %>% enframe() %>% separate(name, into = c("comparison", "field"), sep="\\.") %>%
    spread(key=field, value=value) %>%
    mutate(log2FoldChange = as.numeric(log2FoldChange),
           padj = as.numeric(padj),
           padj_cutoff = !!pthresh,
           absLog2FoldChange_cutoff = !!abslogfcthresh)

  comp = comp %>% mutate(sig = (padj < !!pthresh) & (abs(log2FoldChange) > !!abslogfcthresh)) %>%
    select(comparison, sig, padj, log2FoldChange, everything())

  return(comp)
}

get_PCA <- function(dds, labelled_samples = NULL, color.group, shape.group, variancetransform = T){

  require(DESeq2)
  require(tidyverse)
  require(ggrepel)


  if(variancetransform == T){
    vsd = vst(dds, blind=T)
  } else {
    vsd = dds
  }

  sampleDists <- dist(t(assay(vsd)))
  pcaData <- plotPCA(vsd, intgroup=c(color.group, shape.group), returnData=TRUE) #To do: Make intgroup variables
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  if (is.null(labelled_samples)==T){
    ggplot(pcaData, aes(PC1, PC2, color=get(color.group), shape=get(shape.group))) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() + labs(color = color.group, shape = shape.group)
  } else {
    pcaData = pcaData %>% mutate(labels = if_else(name %in% labelled_samples, name, NULL))
    ggplot(pcaData, aes(PC1, PC2, color=get(color.group), shape=get(shape.group), label=labels)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() + ggrepel::geom_label_repel(size=4, show.legend = F) +
      labs(color = color.group, shape = shape.group)
  }

}

get_PCA_from_matrix <- function(mat, sampledata, labelled_samples = NULL, color.group, shape.group, ntop=500){

  require(tidyverse)
  require(ggrepel)


  rv <- rowVars(mat)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(mat[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)

  if (!all(rownames(sampledata) %in% colnames(mat))) {
    stop("The column names of the input matrix should match the row names of the sample data")
  }

  intgroup = c(color.group, shape.group)
  intgroup.df <- sampledata[, intgroup, drop = FALSE]

  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  } else {
    sampledata[[intgroup]]
  }

  pcaData <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
                        intgroup.df, name = colnames(mat))

  attr(pcaData, "percentVar") <- percentVar[1:2]

  #sampleDists <- dist(t(assay(vsd)))
  #pcaData <- plotPCA(vsd, intgroup=c(color.group, shape.group), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))

  if (is.null(labelled_samples)==T){
    ggplot(pcaData, aes(PC1, PC2, color=get(color.group), shape=get(shape.group))) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() + labs(color = color.group, shape = shape.group)
  } else {
    pcaData = pcaData %>% mutate(labels = if_else(name %in% labelled_samples, name, NULL))
    ggplot(pcaData, aes(PC1, PC2, color=get(color.group), shape=get(shape.group), label=labels)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() + ggrepel::geom_label_repel(size=4, show.legend = F) +
      labs(color = color.group, shape = shape.group)
  }

}

scree_plot = function(res_prcomp, n_pca=20){
  require(tidyverse)
  percentVar <- res_prcomp$sdev^2/sum(res_prcomp$sdev^2) * 100
  names(percentVar) = paste0("PC", seq(1, length(percentVar)))
  percentVar = enframe(percentVar, "PC", "variance")
  percentVar$PC=factor(percentVar$PC, levels=percentVar$PC)
  percentVar %>% slice(1:n_pca) %>% ggplot(aes(x=PC, y=variance)) + geom_bar(stat="identity") +
    ggtitle("Scree plot") + ylab("Percent total variance")

}

gene_report = function(list_reports, ensembl_id, pthresh = 0.05, abslogfcthresh = 0.5){

  require(tidyverse)

  #print(paste("Adjusted p value threshold:", pthresh, ", abs(log2FoldChange) threshold:", abslogfcthresh))

  comp = lapply(list_reports,
                function(x) x[x$ensembl_gene_id==ensembl_id,
                              c("ensembl_gene_id","gene_name","log2FoldChange","padj")]) %>%
    unlist() %>% enframe() %>% separate(name, into = c("comparison", "field"), sep="\\.") %>%
    spread(key=field, value=value) %>%
    mutate(log2FoldChange = as.numeric(log2FoldChange),
           padj = as.numeric(padj),
           padj_cutoff = !!pthresh,
           absLog2FoldChange_cutoff = !!abslogfcthresh)

  comp = comp %>% mutate(sig = (padj < !!pthresh) & (abs(log2FoldChange) > !!abslogfcthresh)) %>%
    select(comparison, sig, padj, log2FoldChange, everything())

  comp = comp %>% arrange(padj)

  return(comp)
}

#### WIP ####

#This swiss army function attempts to aggregate the heat map by duplicate gene names as well as plotting the heatmap
#Appears to have issues in swapping positions for row aggregation colors

complex_heatmap_dedup <- function(vsd, annotated_results, groups_to_plot=levels(vsd$study_group),
                                  row_annot=F, row_scale = FALSE, row_id = "gene_name",
                                  row_size = 8, col_size = 8, dedup_gene_ids = T, title=NULL, ...){

  #This function has some issues with the row annotation getting shuffled if the
  #Deduplication option is used

  require(DESeq2)
  require(ComplexHeatmap)
  require(tidyverse)
  require(scrime)
  require(RColorBrewer)

  stopifnot(row_id %in% c("gene_name", "ensembl_gene_id", "concatenate"))
  stopifnot(groups_to_plot %in% levels(vsd$study_group))
  if (dedup_gene_ids == T & row_id != "gene_name"){
    stop("Deduplication only applicable with gene names")
  }

  genes_to_plot = annotated_results[,c("ensembl_gene_id", "gene_name")]
  samples_to_plot = colnames(vsd)[as.data.frame(colData(vsd)["study_group"])[,1] %in% groups_to_plot]

  vsd = vsd[genes_to_plot$ensembl_gene_id, samples_to_plot]

  mat = assay(vsd)

  if (row_scale==T){
    mat = scrime::rowScales(mat)
  }

  title = if_else(groups_to_plot == levels(vsd$study_group),
                  "All groups",
                  paste(groups_to_plot, collapse="vs"))

  stopifnot(identical(rownames(mat), genes_to_plot$ensembl_gene_id))

  if (row_id == "gene_name"){
    row_labels = genes_to_plot$gene_name
  } else if (row_id == "concatenate"){
    row_labels = paste(genes_to_plot$gene_name, genes_to_plot$ensembl_gene_id, ":")
  } else {
    row_labels = rownames(mat)
  }

  if (dedup_gene_ids == T & sum(duplicated(genes_to_plot$gene_name)) > 0){
    print("Averaging expression for duplicate genes in heatmap matrix")
    mat = rownames_to_column(as.data.frame(mat), "ensembl_gene_id")
    mat = right_join(select(annotated_results, ensembl_gene_id, gene_name), mat, by = "ensembl_gene_id") %>%
      mutate(gene_name = if_else(is.na(gene_name), ensembl_gene_id, gene_name)) %>% select(-ensembl_gene_id)

    mat = summarize_expression_duplicate_ids(mat, "gene_name") #Gene name de-deduplication function
    rownames(mat) = mat[,colnames(mat)=="GeneSymbol"]
    mat = mat[,colnames(mat) != "GeneSymbol"]
    mat = as.matrix(mat)
    row_labels = rownames(mat)
  }

  # Top column annotation
  ann_top = as.data.frame(colData(vsd)["study_group"])

  #Colors must be named lists with named elements
  study_colors=suppressWarnings(colorRampPalette(brewer.pal(length(unique(vsd$study_group)),"Spectral"))(length(unique(vsd$study_group))))
  names(study_colors) = unique(vsd$study_group)

  top_cols=list(study_group=study_colors)

  colTop <- HeatmapAnnotation(df=ann_top, which="col", col = top_cols)

  #Bottom column annotation
  ann_bottom = as.data.frame(colData(vsd)["PAM50"])
  pam_colors = suppressWarnings(colorRampPalette(brewer.pal(length(unique(vsd$PAM50)),"Paired"))(length(unique(vsd$PAM50))))
  names(pam_colors)= unique(vsd$PAM50)
  bottom_cols = list(PAM50 = pam_colors)
  colBottom <- HeatmapAnnotation(df=ann_bottom, which="col", col = bottom_cols)

  # Row annotation
  if (dedup_gene_ids == T){ #Ensure we don't have too many row annotation entries after removing duplicate gene ids
    anno_rows = tibble(ensembl_gene_id = annotated_results$ensembl_gene_id,
                       gene_name = annotated_results$gene_name,
                       gene_type=annotated_results$gene_type) %>%
      filter(!duplicated(gene_name)) %>%
      mutate(Type=case_when(gene_type == "protein_coding" ~ "protein coding",
                            str_detect(string=gene_type,pattern="IG") ~ "IG gene",
                            TRUE ~ "other noncoding")) %>%
      column_to_rownames("ensembl_gene_id") %>% select(Type)
  } else {
    anno_rows = tibble(ensembl_gene_id = annotated_results$ensembl_gene_id, gene_type=annotated_results$gene_type) %>%
      mutate(Type=case_when(gene_type == "protein_coding" ~ "protein coding",
                            str_detect(string=gene_type,pattern="IG") ~ "IG gene",
                            TRUE ~ "other noncoding")) %>%
      column_to_rownames("ensembl_gene_id") %>% select(Type)
  }
  anno_rows$Type = as.factor(anno_rows$Type)
  type_colors = suppressWarnings(colorRampPalette(brewer.pal(length(levels(anno_rows$Type)),"Set1"))(length(levels(anno_rows$Type))))
  names(type_colors) = levels(anno_rows$Type)
  row_colors = list(Type = type_colors)
  rowAnno = HeatmapAnnotation(df=anno_rows, which="row", col=row_colors)

  #Change legend according to whether input is scaled
  if (row_scale==T){
    hlp = list(title="rowscaled vst")
  } else {
    hlp = list(title="vst counts")
  }

  hm = Heatmap(mat, top_annotation = colTop, bottom_annotation = colBottom, left_annotation = rowAnno,
               heatmap_legend_param = hlp, row_labels = row_labels,
               row_names_gp = gpar(fontsize = row_size), column_names_gp = gpar(fontsize = col_size), ...)

}

