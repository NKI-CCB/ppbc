#### Enrichment analysis functions
#Hypergeom tests by Evert Bosdriesz
source(here("src", "enrichment-analysis-functions.R"))


#### Standard report functions ####

annotate_results <- function(results_object, anno_df = gx_annot, immune_df = immune_gene_list,
                             mark_immune = T, mark_outliers = T){
  
  #Takes a res object from deseq, merges it with a gene annotation file, and an optional immune list to create a data frame sorted by adjusted p value
  #Type is an aggregate of immunoglobulin genes from gene type and genes with a known immune function from ImmPort
  require(tidyverse)
  
  anno_res <- results_object %>% as.data.frame() %>% rownames_to_column("ensembl_gene_id") %>%
    right_join(anno_df,., by = "ensembl_gene_id") %>% arrange(padj)
  
  if (mark_immune==T){ #Shows genes that defined as immune based on external list
    anno_res = anno_res %>% mutate(ImmPort_gene = if_else(
      ensembl_gene_id %in% immune_df$ensembl, T, F
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





#' Title Complex Heatmap for Deseq Results
#'
#' @param mat A gene expression matrix with genes as rows and samples as columns. 
#' Rownames must not repeat, and should be found in sig_results$gene_name.
#' @param sampledata A data frame that matches in sample name from colnames(mat) to a phenotype
#' in the heatmap annotation color scheme. Must contain a primary column of interest specified in intgroup
#' and a Type column containing gene biotypes that match the named vector for gene_colors.
#' @param sig_results A data frame created by calling annotate_results() followed by significance() on a dds results object,
#' or the product of shrinkRes(). Only these genes will be plotted. Must contain columns "gene_name" and "Type".
#' @param groups_to_plot A subset of intgroup that will be plotted. Default: levels(sampledata$study_group)
#' @param intgroup The main category of interest and a column in dds
#' @param title A text string to title the heat map
#' @param top_vars A vector of columns from sample data to include in the top heatmap annotation
#' Default =  c("study_group")
#' @param top_colors A named list of of colors to match top_vars.
#' Default: list(study_group=study_colors)
#' @param bottom_vars A vector of columns from sample data to include in the bottom heatmap annotation
#' Default" bottom_vars = c("PAM50")
#' @param bottom_colors A named list of of colors to match top_vars.
#' Default: list(PAM50=bottom_colors)
#' @param row_colors Currently the only valid option is list(Type = gene_colors), where gene_colors is a vector
#' containing a color for each unique value of sig_results$Type.
#' Ex. gene_colors
#' other noncoding immune protein coding        protein coding        immunoglobulin       T cell receptor 
#' "#1B9E77"             "#D95F02"             "#7570B3"             "#E7298A"             "#66A61E" 
#' @param legend_title A string to rename the legend correspoding to the top row annotation. Optional.
#' @param row_scale Whether to mean-center the rows of the heatmap
#' @param maxn_rownames The maximum number of rows to allow before row names are turned off. Default = 50.
#' @param row_size The fontsize of the row (Default: 8)
#' @param col_size The fontsize of the columns (Default: 8)
#' @param show_col_names Whether to show the column names
#' @param ... Additional parameters to pass to Heatmap()
#'
#' @return A heatmap object from ComplexHeatmap
#' @export 
#'
#' @examples deseq_heatmap(mat = dedup_geneEx, sampledata = as.data.frame(colData(variance_stabilized_dds)),
#' sig_results = sig, title = title, groups_to_plot = groups)
#' 
deseq_heatmap = function(mat, sampledata, sig_results,
                         groups_to_plot=levels(sampledata$study_group),
                         title=NULL,
                         intgroup = c("study_group"),
                         top_vars = c("study_group"),
                         top_colors = list(study_group=study_colors),
                         bottom_vars = c("PAM50"),
                         bottom_colors = list(PAM50 = pam_colors),
                         row_colors = list(Type = gene_colors),
                         legend_title = NULL,
                         row_scale = T,
                         maxn_rownames = 50,
                         row_size = 8, col_size = 8,
                         show_col_names = F, ...){
  require(scrime)
  require(ComplexHeatmap)
  require(tidyverse)
  
  stopifnot(groups_to_plot %in% levels(sampledata[,intgroup])) #Now accepts other columns
  stopifnot(identical(colnames(mat), sampledata$sample_name)) #Ensure that a column containing all sample names exists
  
  
  if(is.null(title)==T){
    title = if_else(identical(groups_to_plot,levels(sampledata[, intgroup])),
                    "All groups",
                    paste(groups_to_plot, collapse="vs"))
  }
  
  
  #Reduce genes to significant only
  genestoplot = sig_results$gene_name
  mat = mat[rownames(mat) %in% genestoplot, ]
  
  #Reduce genes and sample data to compared groups only
  #sampledata = sampledata %>% dplyr::filter(study_group %in% groups_to_plot)
  sampledata = sampledata %>% dplyr::filter(!!as.symbol(intgroup) %in% groups_to_plot)
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
  ann_top = sampledata[,top_vars, drop=F]
  
  #Top column annotation
  colTop <- HeatmapAnnotation(df=ann_top, which="col",
                              col = top_colors,
                              annotation_legend_param = list(list(title = legend_title)))
  
  #Bottom column annotation
  ann_bottom = sampledata[,bottom_vars, drop=F]
  colBottom <- HeatmapAnnotation(df=ann_bottom, which="col", col = bottom_colors)
  
  # Row annotation
  anno_rows = sig_results %>%
    select(gene_name, Type) %>%
    distinct() %>%
    filter(!duplicated(gene_name)) %>%
    column_to_rownames("gene_name")
  
  #Essential that the order be the same!
  anno_rows = anno_rows[match(rownames(mat),rownames(anno_rows)), ,drop=F]
  
  rowAnno = HeatmapAnnotation(df=anno_rows, which="row",
                              col=row_colors
                              #col=list(Type = row_colors)
  )
  
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
                             intgroup="study_group", colorby="PAM50", color_fun = scale_color_viridis_d(na.value = "grey50"),
                             title_string = "", path_save_fig = NULL){
  
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
  
  samples_to_keep = colnames(dds)[as.data.frame(colData(dds)[intgroup])[,1] %in% groups_to_plot]
  
  dds.trim = dds[,samples_to_keep]
  
  #Normalizes by size factors
  #print(colnames(colData(dds.trim)))
  #print(intgroup)
  #print(colorby)
  d <- plotCounts(dds.trim, gene=ensembl_gene_id, intgroup=c(intgroup,colorby), normalize=T, returnData=TRUE)
  
  d[,intgroup] = as.factor(d[,intgroup])
  d[,colorby] = as.factor(d[,colorby])
  
  title = paste0(paste0(title_string, ", "),
                 paste(gene_name,ensembl_gene_id, sep=":"), "\npadj: ",
                 formatC(result_df$padj, format = "e", digits = 2),
                 ", log2fc: ", round(result_df$log2FoldChange,2))
  
  dplot = ggplot(d, aes(x=get(intgroup), y=count)) +
    geom_boxplot(alpha=0, show.legend = F, outlier.shape = NA) +
    geom_point(position=position_jitter(w=0.1,h=0), aes(color=get(colorby))) +
    scale_y_log10() + xlab(intgroup) + labs(color=colorby) + ylab("size-factor normalized counts") +
    ggtitle(title) +
    theme_classic() +
    color_fun
  
  if (is.null(path_save_fig)!=T){
    suppressMessages(ggsave(filename=path_save_fig, plot = dplot))
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


#' @title DESeq2 report wrapper
#' 
#' @description Creates a list of important results objects from various DESeq2 analyses for the PPBC project.
#'
#' @param results_object A data frame creating by calling deseq2::results(), deseq2::lfcshrink() or shrinkRes on a deseqdataset.
#' @param dds A deseqdataset with results from DESeq().
#' @param anno_df A data frame containing a gene dictionary with the following columns:
#'  ensembl_gene_id, gene_name, gene_type and description. Additional columns are permitted.
#' @param mark_immune 
#' @param immune_df A data frame containing columns for ensembl ids and gene names of genes with known immune functions.
#' Additional columns are permitted.
#' Default df is from the ImmPort database. See https://www.innatedb.com/redirect.do?go=resourcesGeneLists
#' @param pthresh The adjusted p value (from DESeq2) threshold for significance. Default: 0.05
#' @param absl2fc The absolute log2 fold change (from DESeq2) for significance. Default: 0.5 (equivalent to log fold change of 1.414214)
#' @param variance_stabilized_dds A DESeq2Tranform object created by calling DESeq2::vst(), DESeq2::varianceStabilizingTransformation(),
#' or blind_vst(). This object will be used to create the (row-scaled) heatmap.
#' @param title A string that will be incorporated into the title of plots. Ex: "Involuting vs nulliparous"
#' @param list_gene_sets An optional list of genes sets (for example: GO-BP) in nested format, created by flexsgsea::read_gmt (or equivalent output).
#' Default: gene_sets = list(
#' go_bp = flexgsea::read_gmt(here("data","external","gmt","c5.bp.v7.0.symbols.gmt")),
#' hallmark = flexgsea::read_gmt(here("data","external","gmt","h.all.v7.0.symbols.gmt")),
#' c2_canon = flexgsea::read_gmt(here("data","external","gmt","c2.cp.v7.0.symbols.gmt")),
#' c2_cgp = flexgsea::read_gmt(here("data","external","gmt","c2.cgp.v7.0.symbols.gmt")))
#' Required for fisher's exact pathway analysis, which will be skipped if param is NULL.
#' @param plot_pathway Whether to create a dot plot from Fisher's exact pathway analysis.
#' Will be ignored if list_genes_sets = NULL. Default = T.
#' @param top_vars A vector of columns from colData(variance_stabilized_dds) to include in the top heatmap annotation
#' Default =  c("study_group")
#' @param top_colors A named list of of colors to match top_vars.
#' Default: list(study_group=study_colors)
#' @param bottom_vars A vector of columns from colData(variance_stabilized_dds) to include in the bottom heatmap annotation
#' Default" bottom_vars = c("PAM50")
#' @param bottom_colors A named list of of colors to match bottom_vars.
#' Default: list(PAM50=bottom_colors)
#' @param intgroup The primary column of interest from colData(dds). Will be used to select which groups to include in plots.
#' Default = "study_group"
#' @param colorby A column in colData(dds) which will be colorized in beehive plots.
#' @param groups Which values of \code{dds[,intgroup,]} to include in plots. Can be used to subset intgroup to, for example, plot only
#' the two groups in a pairwise comparison.
#' @param beehive_groups Either "comparison" or "all". A manual toggle to force the groups plotted in beehives to differ from
#' those plotted in the heatmap.
#' @param n_beehive_plots The number of beehive plots to be shown.
#' @param verbose Whether to print progress messages
#' @param ... Additional parameters for deseq_heatmap()
#'
#' @return A list object with the following components:
#' "annotated_results"
#' "sig_threshold"     
#' "significant_genes" 
#' "volcano_plot"      
#' "heatmap"           
#' "beehive_plots"     
#' "pathways"          
#' "pathway_plots"     
#' "cooks_threshold"   
#' "volc_outliers"    
#' "outlier_bees"   
#'
#' @examples rep_inv_nonprbc = deseq_report(ape, dds=dds, absl2fc = 0.5, title="Involution vs nulliparous",
#'                   intgroup="PPBC", groups = c("involuting", "nulliparous"), 
#'                   colorby="survival",
#'                   top_vars = c("PPBC", "survival"),
#'                   top_colors = list(PPBC=ppbc_colors,
#'                                     survival = os_colors),
#'                   bottom_vars = c("ER", "HER2"),
#'                   bottom_colors = list(ER = er_colors,
#'                                        HER2 = her2_colors))

deseq_report = function(results_object, dds, anno_df = gx_annot, mark_immune=T, immune_df=immune_gene_list,
                        pthresh = 0.05, absl2fc = 0.5, variance_stabilized_dds = vsd,
                        title=NULL, list_gene_sets = gene_sets, plot_pathway = T,
                        top_vars = c("study_group"),
                        top_colors = list(study_group=study_colors),
                        bottom_vars = c("PAM50"),
                        bottom_colors = list(PAM50 = pam_colors),
                        intgroup="study_group", colorby="PAM50",
                        groups=levels(colData(dds)[,intgroup]),
                        beehive_groups = "comparison",
                        n_beehive_plots=5, verbose=T, ...){
  require(DESeq2)
  require(tidyverse)
  
  stopifnot(beehive_groups %in% c("all", "comparison"))
  stopifnot(intgroup %in% colnames(colData(dds)))
  
  mega = list()
  
  if(is.null(title==T)){
    title = if_else(identical(sort(groups),
                              sort(levels(colData(dds)[,"study_group"]))), "All groups",
                    paste(groups, collapse=" vs "))
  }
  
  #Merge deseq results object with gene name and type annotation
  if(verbose==T){
    print("Annotating results...")
  }
  
  df = annotate_results(results_object, mark_immune = mark_immune, immune_df = immune_df)
  
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
                     intgroup = intgroup,
                     top_vars = top_vars, bottom_vars = bottom_vars,
                     top_colors = top_colors, bottom_colors = bottom_colors,
                     sig_results = sig, title = title, groups_to_plot = groups, ...)
  
  
  mega$heatmap = hm
  
  #Draw a selection of beehive plots from top hits
  if(verbose==T){
    print(paste("Creating beehive plots from top", n_beehive_plots, "genes..."))
  }
  
  if(beehive_groups == "all"){
    beehive_groups = levels(colData(dds)[,intgroup]) #Plot all study groups in beehive
  } else {
    beehive_groups = groups #Plot only the ones involved in calculating significance
  }
  
  beelist = plot_many_beehives(dds = dds, result_df = head(sig, n_beehive_plots), groups_to_plot = beehive_groups,
                               intgroup=intgroup, colorby=colorby, title_string = title)
  
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
  
  outlier_bees = plot_many_beehives(dds = dds, title_string = title,
                                    result_df = head(outlier_df, 10),
                                    groups_to_plot = beehive_groups,
                                    intgroup=intgroup, colorby=colorby)
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

plot_gene_simple = function(dds, ensembl_gene_id, gene_name, intgroup="study_group", colorby="PAM50", path_save_fig = NULL){
  require(DESeq2)
  require(tidyverse)
  
  #Normalizes by size factors and log transform the data with a pseudocount of 0.5
  d <- plotCounts(dds, gene=ensembl_gene_id, intgroup=c(intgroup,colorby), normalize=T, transform=T, returnData=TRUE)
  
  d[,intgroup] <- as.factor(d[,intgroup])
  d[,colorby] <- as.factor(d[,colorby])
  
  dplot = ggplot(d, aes(x=get(intgroup), y=count, color=get(colorby))) +
    geom_point(position=position_jitter(w=0.1,h=0)) +
    #geom_boxplot(alpha = 0) +
    scale_y_log10() + xlab(intgroup) + labs(color=colorby) +
    ggtitle(paste(gene_name,ensembl_gene_id, sep=":")) +
    theme_minimal()
  
  if (is.null(path_save_fig)!=T){
    suppressMessages(ggsave(filename=path_save_fig, plot = dplot))
  }
  
  return(dplot)
}

#### Comparing two differential expression formulas ----

#' Plot the differences in results from two differential expression formula
#'
#' @param x A results data frame from DESeq2, containing the following fields:
#' `ensembl_gene_id, gene_name, padj, log2FoldChange, cooks_outlier`
#' and a column from which to extract the plot's colors. See: `color.col` 
#' @param y A df as with x, produced using a different differential expression formula
#' @param plot.stat A statistic to plot. Must be either "padj" or "l2fc"
#' @param sig.only Whether to plot only the genes that are significant in at least one comparison
#' @param threshold.padj The adjusted p threshold below which a gene is significant. Default: 0.05
#' @param threshold.l2fc The absolute value of log2 fold change above which a gene is significant. Default: 0.5
#' @param color.col A column from which to extract the plot's colors (default = Type).
#' @param cols A named vector of colors that match the levels of `color.col`
#' @param xlab A label for the x axis. Default: "~batch + PAM50"
#' @param ylab A label for the y axis. Default: "~ER + HER2"
#' @param dp_title A title for the dotplot
#'
#' @return A ggplot dot plot
#'
#' @examples plot_formula_differences(x = bPam$LRT, y = eh$LRT, plot.stat = "padj", sig.only = T)

plot_formula_differences <- function(x, y, plot.stat = "padj", sig.only = F, 
                                     threshold.padj = 0.05, threshold.l2fc = 0.5,
                                     color.col = "Type", cols = colors$gene_colors,
                                     xlab = "~batch + PAM50", ylab ="~ER + HER2", 
                                     dp_title){
  
  require(dplyr)
  require(ggplot2)
  
  stopifnot(plot.stat %in% c("padj", "l2fc"))
  
  
  #Combine x and y
  combo <- dplyr::full_join(dplyr::select(x, ensembl_gene_id, gene_name,
                                          Type, padj, l2fc = log2FoldChange,
                                          outlier = cooks_outlier),
                            dplyr::select(y, ensembl_gene_id, padj,
                                          l2fc = log2FoldChange,
                                          outlier = cooks_outlier),
                            by="ensembl_gene_id") %>%
    dplyr::select(ensembl_gene_id, gene_name, Type,
                  padj.x, padj.y, l2fc.x, l2fc.y, outlier.x, outlier.y)
  
  #Set up a gray color for cases when exceeding the cook's distance causes padj to be set to zero
  combo[,color.col] <- dplyr::case_when(
    combo$outlier.x == T ~ "outlier.x",
    combo$outlier.y == T ~ "outlier.y",
    TRUE ~ combo[,color.col, drop=T]
  )
  
  cols <- c(cols, outlier.x = "gray", outlier.y = "gray")
  
  #Toggleable whether to plot sig genes only
  if(sig.only == T){
    combo <- combo %>%
      dplyr::filter(
        (padj.x < !!threshold.padj & abs(l2fc.x) > !!threshold.l2fc) |
          (padj.y < !!threshold.padj & abs(l2fc.y) > !!threshold.l2fc)
      )
  }
  
  #Store max values of x and y before replacing them with 1
  xmax <- max(combo[,paste0(plot.stat, ".x"), drop=T], na.rm = T)
  ymax <- max(combo[,paste0(plot.stat, ".y"), drop=T], na.rm = T)
  
  #Replace NAs with 1 so that they can be plotted
  combo <- combo %>%
    dplyr::mutate(padj.x = tidyr::replace_na(padj.x, 1),
                  padj.y = tidyr::replace_na(padj.y, 1))
  
  #Set up the title
  dp_title <- paste("Formula comparison for", dp_title)
  if(sig.only==T){
    dp_title = paste0(dp_title, ": ",
                      "significant genes")
  }
  
  
  #Create dotplot
  dotplot <- combo %>%
    ggplot(aes(x = get(paste0(plot.stat, ".x")),
               y = get(paste0(plot.stat, ".y")),
               color = get(color.col))) +
    geom_point() +
    scale_color_manual(values = cols) +
    theme_minimal() + 
    xlab(paste(xlab, plot.stat)) +
    ylab(paste(ylab, plot.stat)) +
    labs(color = color.col) +
    ggtitle(dp_title)
  
  return(dotplot)
}

#' Summarize formula differences
#'
#' @param x A data frame containing all genes, produced by deseq_report
#' @param y A data frame containing all genes, produced by deseq_report
#' @param formula.x The formula used to produce x
#' @param formula.y The formula used to produce y
#' @param threshold.padj The adjusted p value below which a gene is considered significant
#' @param threshold.l2fc The absolute value of log2 fold change above which a gene is considered significant
#' @param colname.x What to call comparison x in plots
#' @param colname.y What to call comparison y in plots
#'
#' @return A list of data frames, containing the following values:
#' `param` A data frame containing metadata about x and y
#' `thresholds` The padj and absl2fc thresholds  
#' `sig.x` A data frame containing all genes significant in x       
#' `sig.y` A data frame containing all genes significant in y              
#' `vennCounts` From limma, a brief overview of the overlap between the data set.
#'    Can be used with limma::vennDiagram
#' `gene.table` A data frame containing an expanded genewise version of vennCounts  
#' `sig.x.not.y` A data frame showing statistics about genes unique to x
#' `sig.y.not.x` A data frame showing statistics about genes unique to y
#' `sig.x.and.y` A data frame showing statistics about genes common to x and y
#'
#' @examples summarize_formula_differences(x = bpam$LRT, y = eh$LRT)
summarize_formula_differences <- function(x, y,
                                          formula.x = "~ batch + PAM50 + study_group",
                                          formula.y = "~ ER + HER2 + study_group",
                                          threshold.padj = 0.05,
                                          threshold.l2fc = 0.5,
                                          colname.x = NULL,
                                          colname.y = NULL) {
  require(tidyverse)
  
  stopifnot(all(sort(unique(x$ensembl_gene_id)) == sort(unique(y$ensembl_gene_id))))
  
  #To do: This function could be a lot more succinct
  res <- list()
  
  #Eventual names of venn diagram for arguments x and y
  #Will be the name of the object passed to x and y if left blank
  if(is.null(colname.x)){
    colname.x = deparse(substitute(x))
  }
  if(is.null(colname.y)){
    colname.y = deparse(substitute(y))
  }
  
  res$param = tibble::tibble(var = c("x", "y"),
                             input = c(colname.x, colname.y),
                             formulas = c(formula.x, formula.y))
  res$thresholds = tibble::tibble(threshold = c("padj", "l2fc"),
                                  value = c(threshold.padj, 
                                            threshold.l2fc))
  x = arrange(x, ensembl_gene_id)
  y = arrange(y, ensembl_gene_id)
  
  #Set up a column to indicate whether a gene in x is significant
  x = x %>% dplyr::mutate(is.sig.x = padj < !!threshold.padj & abs(log2FoldChange) > !!threshold.l2fc)
  x$is.sig.x[is.na(x$is.sig.x)] <- F
  gene.x = x %>% filter(is.sig.x == T) %>% pull(ensembl_gene_id)
  
  #Set up a column to indicate whether a gene in y is significant  
  y = y %>% dplyr::mutate(is.sig.y = padj < !!threshold.padj & abs(log2FoldChange) > !!threshold.l2fc)
  y$is.sig.y[is.na(y$is.sig.y)] <- F
  gene.y = y %>% filter(is.sig.y == T) %>% pull(ensembl_gene_id)
  
  #Base diffex lists
  res$sig.x = x %>% dplyr::filter(is.sig.x == T) %>%
    select(gene_name, padj, log2FoldChange, everything())
  res$sig.y = y %>% dplyr::filter(is.sig.y == T) %>%
    select(gene_name, padj, log2FoldChange, everything())
  
  # Count unique and common hits with limma::Venncounts 
  # Requires a data frame that has
  cb <- data.frame(row.names = sort(x$ensembl_gene_id),
                   x = x %>% arrange(ensembl_gene_id) %>% pull(is.sig.x),
                   y = y %>% arrange(ensembl_gene_id) %>% pull(is.sig.y)
  )
  #In these times of social distancing we really can't have too many sanity checks
  stopifnot(nrow(res$sig.x) == sum(cb$x) & nrow(res$sig.y) == sum(cb$y))
  #return(cb)
  
  # Replace columns with the symbols passed to x and y
  v = cb
  colnames(v) = c(colname.x, colname.y)
  v = limma::vennCounts(v)
  res$vennCounts <- v
  
  #Can't easily be saved, so just call the function later
  #res$venn <- limma::vennDiagram(v)
  
  #More logicals for easy retrieval
  cb = cb %>% rownames_to_column("ensembl_gene_id") %>%
    mutate(x.not.y = x == T & y == F,
           y.not.x = y == T & x == F,
           x.and.y = x == T & y == T)
  #  return(cb)
  stopifnot(all(sum(rowSums(cb[,4:6]) == 1)))
  
  cb = right_join(select(x, ensembl_gene_id, gene_name, type = Type), cb,
                  by = "ensembl_gene_id")
  
  #Significance overlap table
  res$gene.table = cb
  #return(cb)
  
  #Retrieve information for genes unique to x, including their results in y
  res$sig.x.not.y = full_join(
    select(filter(x, ensembl_gene_id %in% cb$ensembl_gene_id[cb$x.not.y]),
           padj.x = padj, l2fc.x = log2FoldChange, type = Type, everything()),
    select(filter(y, ensembl_gene_id %in% cb$ensembl_gene_id[cb$x.not.y]),
           padj.y = padj, l2fc.y = log2FoldChange, ensembl_gene_id, is.sig.y),
    by="ensembl_gene_id") %>%
    mutate(padj.diff = padj.x - padj.y,
           l2fc.diff = l2fc.x - l2fc.y) %>%
    select(gene_name, padj.x, padj.y, padj.diff, l2fc.x, l2fc.y, l2fc.diff,  is.sig.x, is.sig.y, everything())
  stopifnot(all(res$sig.x.not.y$is.sig.x == T) & all(res$sig.x.not.y$is.sig.y == F))
  
  #Retrieve information for genes unique to y, including their results in x
  res$sig.y.not.x = full_join(
    select(filter(x, ensembl_gene_id %in% cb$ensembl_gene_id[cb$y.not.x]),
           padj.x = padj, l2fc.x = log2FoldChange, type = Type, everything()),
    select(filter(y, ensembl_gene_id %in% cb$ensembl_gene_id[cb$y.not.x]),
           padj.y = padj, l2fc.y = log2FoldChange, ensembl_gene_id, is.sig.y),
    by="ensembl_gene_id") %>%
    mutate(padj.diff = padj.y - padj.x,
           l2fc.diff = l2fc.y - l2fc.x) %>%
    select(gene_name, padj.x, padj.y, padj.diff, l2fc.x, l2fc.y, l2fc.diff, is.sig.x, is.sig.y,  everything())  
  stopifnot(all(res$sig.y.not.x$is.sig.x == F) & all(res$sig.y.not.x$is.sig.y == T))
  
  #Retrieve information for genes significant in both datasets
  res$sig.x.and.y = full_join(
    select(filter(x, ensembl_gene_id %in% cb$ensembl_gene_id[cb$x.and.y]),
           padj.x = padj, l2fc.x = log2FoldChange, type = Type, everything()),
    select(filter(y, ensembl_gene_id %in% cb$ensembl_gene_id[cb$x.and.y]),
           padj.y = padj, l2fc.y = log2FoldChange, ensembl_gene_id, is.sig.y),
    by="ensembl_gene_id") %>%
    mutate(padj.diff = padj.x - padj.y,
           l2fc.diff = l2fc.x - l2fc.y) %>%
    select(gene_name, padj.x, padj.y, padj.diff, l2fc.x, l2fc.y, l2fc.diff, is.sig.x, is.sig.y, everything())
  stopifnot(all(res$sig.x.and.y$is.sig.x == T) & all(res$sig.x.and.y$is.sig.y == T))
  
  #Reporting
  print(paste("Forumla x:", formula.x))
  print(paste("Forumla y:", formula.y))
  print(paste0("Number of sig ensembl IDs with formula x: ",
               nrow(res$sig.x)))
  #print(paste("IDs unique to x:", v[, "Counts"][3]))
  print(paste("IDs unique to x:", length(unique(res$sig.x.not.y$ensembl_gene_id))))
  print(paste0("Number of sig ensembl IDs with formula y: ",
               nrow(res$sig.y)))
  #print(paste("IDs unique to y:", v[,"Counts"][4]))
  print(paste("IDs unique to y:", length(unique(res$sig.y.not.x$ensembl_gene_id))))
  
  return(res)
}

#' Dotplot of formula differences
#'
#' @description Similar to plot_formula_differences, but designed to work with the results object
#' of summarize_formula_differences
#'
#' @param report_df A data frame in the list produced but summarize_formula_differences
#' @param stat Either "padj" or "l2fc"
#' @param color The column corresponding to the color
#' @param color_vec A named vector of colors where the names match the levels of the color column
#' @param padj.x.label.thresh Show gene names if padj.x above this threshold
#' @param padj.y.label.thresh Show gene names if padj.y above this threshold
#' @param show.inverted.signs Whether to show names of genes that switch signs in l2fc plots
#'
#' @return A ggplot
#'
#' @examples formula_difference_dotplots(lrt_diff$sig.x.not.y, stat = "padj", padj.y.label.thresh = 0.25)
formula_difference_dotplots <- function(report_df, stat = c("padj", "l2fc"),
                                        color = "type",
                                        color_vec = colors$gene_color,
                                        padj.x.label.thresh = Inf,
                                        padj.y.label.thresh = Inf,
                                        show.inverted.signs = F){
  
  #require(tidyverse)
  
  stat <- match.arg(stat)
  rep.name <- deparse(substitute(report_df))
  
  stat.x <- paste0(stat, ".x")
  stat.y <- paste0(stat, ".y")
  
  p <- report_df %>% 
    ggplot(aes(x = get(stat.x),
               y = get(stat.y),
               color = get(color))) +
    geom_point() +
    scale_color_manual(values = color_vec) +
    theme_minimal() +
    ggtitle(paste("Adjusted p of genes significant for", rep.name)) +
    xlab(stat.x) + ylab(stat.y) +
    labs(color = color)
  
  if(stat == "padj"){
    p + geom_hline(yintercept = 0.05, linetype = "dashed") +
      geom_vline(xintercept = 0.05, linetype = "dashed") +
      ggrepel::geom_label_repel(aes(label=gene_name),
                                show.legend = F, 
                                data = filter(report_df,
                                              ((padj.y > !!padj.y.label.thresh) | (padj.x > !!padj.x.label.thresh))
                                ))
  } else {
    p <- p + 
      geom_hline(yintercept = 0, color="gray") +
      geom_vline(xintercept = 0, color="gray")
    if (show.inverted.signs == T){
      p + ggrepel::geom_label_repel(aes(label=gene_name),
                                    show.legend = F,
                                    data = report_df %>%
                                      filter((l2fc.y > 0 & l2fc.x < 0) |
                                               (l2fc.y < 0 & l2fc.x > 0)))
    } else {
      p
    }
    
  }
  
}

#' Deviance from Deseq Datasets
#'
#' @description Unbeknowst to us, DESeq actually does report some genewise deviance measures.
#' However, the statistic reported is -2 times the log likelihood and is intended for calculating LRTs rather than as a goodness of fit metric.
#' For the usual interpretation of deviance, the deviances should be transformed via this funciton.
#' [Source](https://support.bioconductor.org/p/117448/)
#'
#' @param dds A DESeq dataset, after DESeq() has been called.
#' @param print_median Whether to print the median
#'
#' @return A named vector of genewise deviances.
#'
#' @examples dds_dev(dds.x)
dds_dev <- function(dds, print_median = T){
  dev = mcols(dds)$deviance - 2*rowSums(dnbinom(counts(dds), mu=counts(dds), size=1/dispersions(dds), log=TRUE))
  if(print_median){
    print(paste("Median:", round(median(dev),2)))
  }
  return(dev)
}


#' Show mean deviance trend for DDS
#'
#' @param dds A dds
#' @param annot A data frame containing both gene names and ensembl gene ids
#' @param title The title of the dds for the plot
#' @param label.dev.thresh Genes above this and basemean will be labelled
#' @param label.basemean.thresh Genes above this and dev.thresh will be labelled 
#'
#' @return A ggplot dotplot of mean vs deviance
#' @examples  dds_mean_deviance_trend(dds)
dds_mean_deviance_trend <- function(dds, annot = gx_annot,
                                    title = deparse(substitute(dds)),
                                    label.dev.thresh = Inf,
                                    label.basemean.thresh = Inf){
  
  dev = dds_dev(dds, print_median = F)
  base_mean = mcols(dds)$baseMean
  
  mdev_title = paste("Mean - deviance trend for", deparse(substitute(dds)))
  
  df = tibble(ensembl_gene_id = names(dev),
              deviance = dev,
              baseMean = base_mean) %>%
    right_join(select(annot, ensembl_gene_id, gene_name),
               ., by = "ensembl_gene_id")
  
  mdev <- df %>%
    ggplot(aes(x = baseMean, y = deviance, label = gene_name)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    theme_minimal() +
    ggrepel::geom_label_repel(
      data = filter(df,
                    deviance > !!label.dev.thresh & baseMean > !!label.basemean.thresh)
    ) +
    ggtitle(mdev_title)
  
  mdev
}

#' Compare deviance histograms for two dds
#'
#' @param dds.x A dds
#' @param dds.y A different dds
#' @param name.x The name of dds.x for plotting
#' @param name.y The name of dds.y for plotting
#'
#' @return A ggplot histogram with the deviances of x and y
#' @examples compare_deviance_histogram(dds.x, dds.y)
compare_deviance_histogram <- function(dds.x, dds.y,
                                       name.x = deparse(substitute(dds.x)),
                                       name.y = deparse(substitute(dds.y))){
  
  dev.x <- dds_dev(dds.x, print_median = F)
  dev.y <- dds_dev(dds.y, print_median = F)
  
  df <- rbind(tibble(ensembl_gene_id = names(dev.x),
                     deviance = dev.x,
                     name = name.x),
              tibble(ensembl_gene_id = names(dev.y),
                     deviance = dev.y,
                     name = name.y)
  )  
  
  df %>%
    ggplot(aes(x = deviance, fill = name)) +
    geom_histogram(alpha = 0.3, binwidth = 100, position = "identity") +
    theme_minimal()+
    ggthemes::scale_fill_colorblind()+
    ggtitle(paste("Histogram of deviances for\n", name.x,
                  paste0("(median: ", round(median(dev.x),2), ")"),
                  "\n", name.y,  paste0("(median: ", round(median(dev.y),2), ")")))
  
  
}

#' Compare mean deviance trends from 2 dds
#'
#' @param dds.x A dds
#' @param dds.y A different dds
#' @param annot A data frame containing both gene names and ensembl gene ids
#' @param name.x Name dds.x for plotting
#' @param name.y Name dds.y for plotting
#' @param label.dev.thresh Genes above this and basemean will be labelled
#' @param label.basemean.thresh Genes above this and dev.thresh will be labelled 
#'
#' @return A gg dot plot with mean on the x axis and variance on the y axis,
#' with two smoothed lines representing dds.x and dds.y
#'
#' @examples compare_mean_deviance_trend(dds.x, dds.y)
compare_mean_deviance_trend <- function(dds.x, dds.y, annot = gx_annot,
                                        name.x = deparse(substitute(dds.x)),
                                        name.y = deparse(substitute(dds.y)),
                                        label.dev.thresh = Inf,
                                        label.basemean.thresh = Inf){
  
  dev.x = dds_dev(dds.x, print_median = F)
  dev.y = dds_dev(dds.y, print_median = F)
  base_mean.x = mcols(dds.x)$baseMean
  base_mean.y = mcols(dds.y)$baseMean
  
  
  df = rbind(tibble(ensembl_gene_id = names(dev.x),
                    deviance = dev.x,
                    baseMean = base_mean.x,
                    name = name.x),
             tibble(ensembl_gene_id = names(dev.y),
                    deviance = dev.y,
                    baseMean = base_mean.y,
                    name = name.y)) %>%
    right_join(select(annot, ensembl_gene_id, gene_name),
               ., by = "ensembl_gene_id")
  
  df %>%
    ggplot(aes(x = baseMean, y = deviance,label = gene_name)) +
    geom_point(alpha = 0.5, aes(shape = name), color = "gray") +
    geom_smooth(aes(color = name)) +
    scale_y_log10() +
    scale_x_log10() +
    theme_minimal() +
    ggrepel::geom_label_repel(show.legend = F,
                              data = filter(df,
                                            deviance > !!label.dev.thresh & baseMean > !!label.basemean.thresh)
    ) +
    ggtitle(paste("Mean-deviance trends for", name.x, "and", name.y)) +
    ggthemes::scale_color_few()
  
  
}

plot_high_deviance_gene <- function(dds, gene_name, ensembl_gene_id, dev_df = devs.both,
                                    intgroup.x = "batch", colorby.x = "PAM50",
                                    intgroup.y = "ER", colorby.y = "HER2",
                                    formula.x = "formula x", formula.y = "formula y"){
  
  dev.diff = dev_df %>% filter(ensembl_gene_id == !!ensembl_gene_id) %>% pull(dev.diff)
  
  ggpubr::ggarrange(
    plot_gene_simple(dds, ensembl_gene_id = ensembl_gene_id, gene_name = gene_name, intgroup=intgroup.x, colorby=colorby.x)+
      ggtitle(formula.x) + ylab("size factor normalized counts"),
    plot_gene_simple(dds, ensembl_gene_id = ensembl_gene_id, gene_name = gene_name, intgroup=intgroup.y, colorby=colorby.y) +
      ggtitle("formula y") + ylab("size factor normalized counts"),
    ncol =2) %>%
    ggpubr::annotate_figure(top = paste("Gene", gene_name, "dev y - x =", signif(dev.diff,2)))
}

#### Currently unused ####

plot_cooks <- function(dds, colorby="study_group", method="box", outlier.shape=NA, groupwise_threshold=F, returnOutliers=F,...){
  require(tidyverse)
  require(ggrepel)
  require(RColorBrewer)
  
  stopifnot(method %in% c("box", "upperquartile"))
  if (groupwise_threshold == T & method != "upperquartile"){
    stop("Use groupwise threshold with method = upperquartile.")
  }
  if (returnOutliers == T & groupwise_threshold == F){
    stop("Return outliers only works when groupwise_threshold is TRUE.")
  }
  
  sampledata <- data.frame(sample = rownames(as.data.frame(colData(dds))),
                           group = as.data.frame(colData(dds))[,colorby],
                           stringsAsFactors = F)
  
  #Get the Cooks distances
  cooks <- log10(assays(dds)[["cooks"]])
  
  #Melt into tidy data frame
  cooks_df <- cooks %>% as.data.frame() %>% rownames_to_column("gene_id") %>% gather(key=sample,-gene_id,value="cooks")
  
  #Add study group info
  cooks_df <- cooks_df %>% left_join(., sampledata, by="sample")
  
  #Sort by the color variable so all the colors are plotted together
  cooks_df <- cooks_df %>% arrange(group)
  
  
  #Create a box plot of all cooks distances by samples, including or excluding the outliers depending on outlier.shape
  if (method=="box"){
    graph = cooks_df %>%
      ggplot(aes(x=factor(sample, levels=unique(sample)), y=cooks, color=group)) +
      geom_boxplot(outlier.shape = outlier.shape) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(color=colorby) + xlab("sample") +
      scale_colour_brewer(type = "qual", palette = "Set1")
  }
  
  
  #Plot just the upper quartile of the cooks distances as a single point
  #Allows more readable labels
  if (method=="upperquartile"){
    cooksQuartile = cooks_df %>% group_by(sample) %>% summarise(upper_cook = quantile(cooks)[4]) #Get the upper quartile
    cooksQuartile = cooksQuartile %>% left_join(.,sampledata, by="sample")
    
    #Sort by group to have colors plotted together
    cooksQuartile = cooksQuartile %>% arrange(group)
    
    if (groupwise_threshold==T){
      
      #Set a threshold of two standard deviations above the mean upper quartile for Cook's distance
      cooksThreshold = cooksQuartile %>% group_by(group) %>%
        summarise(meanUC=mean(upper_cook), stdevUC=sd(upper_cook, na.rm=T)) %>% mutate(threshold=meanUC+2*stdevUC)
      
      cooksThreshold = left_join(cooksQuartile,cooksThreshold, by="group")  %>%
        mutate(label = if_else(upper_cook > threshold, sample, NULL)) #Label those samples above the threshold
      
      #First plot the UQ cooks distance by color group
      graph = cooksThreshold %>% ggplot(aes(x=factor(sample, level=unique(sample)), y = upper_cook, color=group, label=label)) +
        geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(color=colorby, title="Labelled samples exceed groupwise Cooks threshold of 2sd + mean(UQ)") + xlab("sample") +
        scale_colour_brewer(type = "qual", palette = "Set1")
      
      #Add labels for samples above the threshold
      graph = graph + ggrepel::geom_label_repel(size=4, show.legend = F)
      
      #store the outliers as a vector
      outliers = cooksThreshold %>% filter(!is.na(label)) %>% select(sample, group)
      
    } else {
      
      #Produce a graph of upper quartiles cooks distance with no samples labeled
      graph = cooksQuartile %>% ggplot(aes(x=factor(sample, level=unique(sample)), y = upper_cook,color=group)) +
        geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        labs(color=colorby) + xlab("sample") +
        scale_colour_brewer(type = "qual", palette = "Set1")
    }
  }
  if (returnOutliers == T){
    return(outliers)
  } else {
    return(graph)
  }
}

voom_from_deseq <- function(dds){
  
  require(DESeq2)
  require(edgeR)
  require(limma)
  
  rawcount = counts(dds, normalized=F)
  dge <- edgeR::DGEList(rawcount)
  dge <- edgeR::calcNormFactors(dge, method='TMM')
  voom <- limma::voom(dge, plot=F)
  vt <- t(voom$E) #Must be transposed for genefu compatibility
  return(vt)
}



