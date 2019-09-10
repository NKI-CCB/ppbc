
annotate_results <- function(results_object, anno_df = gx_annot, mark_immune = T, immune_list = immune_gene_list){
  require(tidyverse)

  anno_res <- results_object %>% as.data.frame() %>% rownames_to_column("ensembl_gene_id") %>%
    right_join(anno_df,., by = "ensembl_gene_id") %>% arrange(padj)

  if (mark_immune==T){
    anno_res = anno_res %>% mutate(immune_gene = if_else(
      ensembl_gene_id %in% immune_list$ensembl |
        str_detect(string = description, "immuno") |
        str_detect(string = description, "immune") |
        str_detect(string = description, "interleukin"),
      TRUE, FALSE
    )) %>% #When description is NA, immune_gene is also NA
      mutate(immune_gene = replace_na(immune_gene, FALSE))
  }

  return(anno_res)
}


significance <- function(annotated_results, pthresh = 0.05, absl2fc = 1){
  require(tidyverse)
  sig = annotated_results %>% filter(padj < !!pthresh & (abs(log2FoldChange) > !!absl2fc))
  return(sig)
}

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

volcano_plot <- function(annotated_results, title, path_save_fig = NULL, pthresh = 0.05, absl2fc=1, shape_col=NULL, dedup_ids = T){

  require(tidyverse)
  require(ggrepel)
  require(here)

  volcdf <- annotated_results %>% mutate(NegLog10FDR = -log10(padj),
                                         Significant = if_else(
                                           (padj < !!pthresh & abs(log2FoldChange) > !!absl2fc),
                                           TRUE, FALSE),
                                         Type=case_when(
                                           gene_type == "protein_coding" ~ "protein coding",
                                           str_detect(string=gene_type,pattern="IG") ~ "IG gene",
                                           TRUE ~ "other noncoding"),
                                         Color = if_else(
                                           Significant == T, Type, "n.s."
                                         )) %>%
    filter(!is.na(Significant))

  if (dedup_ids == T){
    print("Summarizing -log10FDR and log2FoldChange by median for duplicate gene names in volcano plot")
    volcdf = volcdf %>% group_by(gene_name, gene_type, immune_gene, Significant, Type, Color) %>%
      summarise(log2FoldChange = median(log2FoldChange), NegLog10FDR = median(NegLog10FDR)) %>%
      arrange(desc(NegLog10FDR))
  }

  colors = c(`IG gene` = "red", `n.s.` = "gray",`other noncoding`="purple",`protein coding`="darkgreen")

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



complex_heatmap <- function(vsd, annotated_results, groups_to_plot=levels(vsd$study_group),
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


plot_gene_beehive = function(dds, result_df, groups_to_plot = levels(colData(dds)[,intgroup]),
                             intgroup="study_group", colorby="PAM50", path_save_fig = NULL){

  require(DESeq2)
  require(tidyverse)

  if(nrow(result_df) > 1){
    warning("More than one result provided, plotting first row")
    result_df = head(result_df, 1)
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
                        pthresh = 0.05, absl2fc = 1, variance_stabilized_dds = vsd,
                        intgroup="study_group", colorby="PAM50", shape_volc_col = "immune_gene",
                        groups=levels(colData(dds)[,intgroup]),
                        n_beehive_plots=5, verbose=F, ...){
  require(DESeq2)
  require(tidyverse)

  mega = list()

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

  title = if_else(identical(sort(groups), sort(levels(colData(dds)[,"study_group"]))), "All groups - Interpret with caution",
                  paste(groups, collapse=" vs "))

  vp = volcano_plot(annotated_results = df, title=title, pthresh = pthresh, absl2fc = absl2fc, shape_col = shape_volc_col)

  mega$volcano_plot = vp

  #Draw heatmap
  if(verbose==T){
    print("Creating heatmap...")
  }

  #Show row names if fewer than 50 genes plotted
  rn = nrow(sig) <= 50

  hm = complex_heatmap(vsd = variance_stabilized_dds, annotated_results = sig, groups_to_plot=groups, row_scale = T,
                       show_column_names=F, show_row_names=rn, title = title)

  mega$heatmap = hm

  #Draw a selection of beehive plots from top hits
  if(verbose==T){
    print(paste("Creating beehive plots from top", n_beehive_plots, "genes..."))
  }
  beelist = plot_many_beehives(dds = dds, result_df = head(sig, n_beehive_plots), groups_to_plot = groups)

  mega$beehive_plots = beelist

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
