####02_QC_salmon####

retrieve_fastqc_module = function(qc, sample_annot, module.type, dir=here("results/fastqc"), failorwarn = c("FAIL", "WARN"),
                                  min.lib.size = 10^6.7){
  require(tidyverse)
  require(fastqcr)

  qc = qc %>% filter(module == !!module.type) %>% filter(status %in% failorwarn) %>%
    filter(as.integer(tot.seq) > min.lib.size)

  qc = left_join(qc, select(sample_annot, sample_id, sample_name), by= c("sample"="sample_id"))

  file.paths = file.path(dir, qc$sample)

  or = fastqcr::qc_read_collection(file = file.paths,
                                   sample_names = qc$sample_name,
                                   modules = module.type,
                                   verbose=F)[[1]]
  or = left_join(or, sample_annot, by = c("sample"="sample_name"))

  or = left_join(or, select(qc, -module), by=c("sample_id"="sample"))

  return(or)
}


add_qc_to_metadata = function(metadata, qcdata, module.type){
  qcdata = qcdata %>% filter(module == !!module.type) %>% rename(status.module=status)
  newdata = left_join(metadata, select(qcdata, sample, status.module), by=c("sample_id"="sample"))
  colnames(newdata)[length(colnames(newdata))] = str_to_lower(str_replace_all(module.type, " ", "_"))
  return(newdata)
}

process_blast_results <- function(blast_results, annot){
  blast_results = blast_results %>% select(id = X1, refseq_id = X2)
  results = right_join(refseq_db,blast_results, by="refseq_id")
  return(results)
}

get_correlations_for_patient_samples <- function(pref){
  stopifnot(pref %in% c(sample_annot$patient_ref, "random"))
  if (pref %in% sample_annot$patient_ref) {
    cors <- cor(
      tx$counts[, filter(sample_annot, patient_ref %in% pref)$sample_id],
      method = "spearman"
    )
  } else{
    samples <- sample(sample_annot$sample_id, 2)
    cors <- cor(tx$counts[, samples], method = "spearman")
  }
  cors[lower.tri(cors, diag = T)] <- NA
  cors %>%
    as_data_frame(rownames = "sample_1") %>%
    gather(sample_2, spearman, -sample_1) %>%
    filter(!is.na(spearman)) %>%
    mutate(patient_ref = pref) %>%
    select(patient_ref, sample_1, sample_2, spearman)
}

subset_tximport <- function(txi, cols){
  lapply(txi, function(x) if ( is.matrix(x) ) return(x[, cols]) else return(x))
}

####Differential expression####

annotate_results <- function(results_object){
  anno_res <- results_object %>% as.data.frame() %>% rownames_to_column() %>%
    right_join(gx_annot,., by = c("gene_id"="rowname")) %>%
    arrange(padj)
  return(anno_res)
}

####Visualization functions####
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

get_PCA <- function(dds, labelled_samples = NULL, color.group, shape.group){

  require(DESeq2)
  require(tidyverse)
  require(ggrepel)

  vsd = vst(dds, blind=FALSE)

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


complex_heatmap <- function(vsd, annotated_results, groups_to_plot=levels(vsd$study_group),
                            row_annot=F, row_scale = FALSE, row_id = "gene_name",
                            row_size = 8, col_size = 8, dedup_gene_ids = T, title=NULL, ...){

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

    mat = summarize_expression_duplicate_ids(mat, "gene_name")
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


  if (row_scale==T){
    hlp = list(title="rowscaled vst")
  } else {
    hlp = list(title="vst counts")
  }

  hm = Heatmap(mat, top_annotation = colTop, bottom_annotation = colBottom, left_annotation = rowAnno,
               heatmap_legend_param = hlp, row_labels = row_labels,
               row_names_gp = gpar(fontsize = row_size), column_names_gp = gpar(fontsize = col_size), ...)
  #hm = draw(hm, column_title = title)
  #return(hm)
}



####03 IHC and PAM50####
plot_receptor = function(data, receptor, x.column="molecular_subtype", y.column="normCount"){

  require(tidyverse)
  require(ggbeeswarm)

  data %>% filter(gene_name == !!receptor) %>%
    ggplot(aes(x=get(x.column), y = get(y.column), color = get(receptor))) +
    ggbeeswarm::geom_beeswarm() +
    xlab(x.column) + ylab(y.column) + labs(color=paste(receptor, "IHC status")) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(paste(receptor, "variance stabilized transformed mRNA expression"))
}

deduplicate_annotation <- function(annot){

  require(genefu)
  data('pam50', package='genefu', envir=environment())
  pam50.names <- rownames(pam50$centroids)
  dups <- annot[which(duplicated(annot$probe)), "Gene.symbol"]
  print("PAM50 genes with multiple ensembl ids:")
  print(table(pam50.names %in% dups))

  annot <- annot[order(annot$EntrezGene.ID),]
  annot <- annot[!(duplicated(annot$probe)),]

  #Complains about being a tibble
  annot <- as.data.frame(annot)
  rownames(annot) <- annot$probe

  stopifnot(sum(duplicated(annot$probe))==0)

  return(annot)
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

barplot_from_table = function(table, xlab, ylab){
  df = as.data.frame(table)
  df %>% ggplot(aes(x=Var1, y = Freq, fill=Var2)) +
    geom_bar(stat="identity") + xlab(xlab) + labs(fill=ylab)
}

get_max_prob = function(prob_matrix, samples=rownames(prob_matrix)){
  mat = prob_matrix[rownames(prob_matrix) %in% samples, ]
  max_prob = apply(mat, 1, function(x) max(x, na.rm=T))

  colmax = colnames(mat)[max.col(mat)]

  df = tibble(sample_name = names(max_prob),
              prediction = colmax,
              probability = max_prob)

  return(df)
}

####04_tumor_purity####

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

#### Batch effects ####

#Using a variance stabilizing transformation that is fully blind to experimental design.

blind_vst = function(dds){
  design(dds) = formula(~ 1)
  vsd = vst(dds, blind=T)
  mat = assay(vsd)
  return(mat)
}


#### Cooks distance ####

# Examine Cooks distance outliers


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

  #Spare ourselves the nightmare of trying to use dplyr with variable column names
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


#### Differential expression ####


plot_gene_simple = function(dds, ensembl_gene_id, gene_name, intgroup="study_group", colorby="PAM50", path_save_fig = NULL){
  require(DESeq2)
  require(tidyverse)

  #Normalizes by size factors and log transform the data with a pseudocount of 0.5
  d <- plotCounts(dds, gene=ensembl_gene_id, intgroup=c(intgroup,colorby), normalize=T, transform=T, returnData=TRUE)

  dplot = ggplot(d, aes(x=get(intgroup), y=count, color=get(colorby))) +
    geom_point(position=position_jitter(w=0.1,h=0)) +
    scale_y_log10() + xlab(intgroup) + labs(color=colorby) +
    ggtitle(paste(gene_name,ensembl_gene_id, sep=":"))

  if (is.null(path_save_fig)!=T){
    suppressMessages(ggsave(filename=path_save_fig, plot = volcplot))
  }

  return(dplot)
}


plot_gene_beehive = function(dds, result_df, groups_to_plot = NULL,
                             intgroup="study_group", colorby="PAM50", path_save_fig = NULL){

  require(DESeq2)
  require(tidyverse)

  if(nrow(result_df) > 1){
    warning("More than one result provided, plotting first row")
    result_df = head(result_df, 1)
  }

  ensembl_gene_id = result_df$ensembl_gene_id
  gene_name = result_df$gene_name

  if (is.null(groups_to_plot)==T){
    groups_to_plot = levels(colData(dds)[,intgroup])
  }

  comparison = paste(groups_to_plot, collapse = " vs ")

  samples_to_keep = colnames(dds)[as.data.frame(colData(dds)[intgroup])[,1] %in% groups_to_plot]

  dds.trim = dds[,samples_to_keep]

  #Normalizes by size factors
  d <- plotCounts(dds.trim, gene=ensembl_gene_id, intgroup=c(intgroup,colorby), normalize=T, returnData=TRUE)

  title = paste0(paste(gene_name,ensembl_gene_id, sep=":"), ", ", comparison, " padj: ", formatC(result_df$padj, format = "e", digits = 2),
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