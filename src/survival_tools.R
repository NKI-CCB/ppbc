#' Retrieve ntiles of a given gene from expression data
#'
#' @param gene_id An ensembl id or gene symbol (name)
#' @param id_type Either "ensembl" or "symbol"
#' @param gene_dict A data frame containing a column "gene name" with gene symbols and a column "ensembl_gene_id" with ensembl gene ids
#' @param geneEx A gene expression matrix with ensembl IDs as row names. Should be tmm/log2 normalized for compatibility with other functions (this function will work with any metric).
#' @param n The number of desired ntiles. Default = 3.
#' @param labels Text labels for ntiles. Default = c("low", "medium", "high")
#'
#' @return A data frame showing the ntiles for the gene in question
#' @export
#'
#' @examples gene_ntile("ENSG00000000003", id_type = "ensembl"), gene_ntile("MUC2", id_type = "symbol")

gene_ntile <- function(gene_id, id_type, gene_dict = gx_annot, geneEx = ens_mat,
                       n = 3, labels = c("low", "medium", "high")){
  
  stopifnot(id_type %in% c("ensembl", "symbol"))
  stopifnot(length(labels) == n)
  
  if(id_type == "symbol"){
    id <- gene_dict[gene_dict$gene_name == gene_id,]$ensembl_gene_id
    name <- gene_id
    #return(gene_id)
  }
  
  if(id_type == "ensembl"){
    id <- gene_id
    name <- gene_dict[gene_dict$ensembl_gene_id == gene_id,]$gene_name
  }
  
  gene <- ens_mat[rownames(ens_mat) %in% id, , drop = F]
  if(nrow(gene) == 0){
    stop("Gene not found in provided count matrix")
  }
  
  if(nrow(gene) > 1){
    warning("Multiple ensembl ids mapped to single gene symbol. Averaging expression.")
    gene <- t(as.matrix(colMeans(gene)))
    #return(gene)
  }
  
  rownames(gene) <- name #Not strictly necessary, useful for debugging
  
  if(length(id) > 1){
    id = paste(id, collapse = ";")
  }
  
  gene_df <- tibble::tibble(sample_name = colnames(gene),
                            gene_symbol = name,
                            ensembl_gene_id = id,
                            tmm = gene[1,],
                            ntile = dplyr::ntile(gene[1,], n),
                            labels = factor(dplyr::ntile(gene[1,], n), labels = labels)
  )
  
  #gene_df <- dplyr::mutate(gene_df, labels = as.character(labels)) #Or survminer will complain
  
  return(gene_df)
}



#' Plot Kaplan-Meier curves of ntiles of a given gene
#'
#' @param gene_id An ensembl id or gene symbol
#' @param id_type Must be either "symbol" or "ensembl"
#' @param gene_dict A data frame containing a column "gene name" with gene symbols and a column "ensembl_gene_id" with ensembl gene ids
#' @param geneEx A gene expression matrix with ensembl IDs as row names. Should be tmm/log2 normalized for compatibility with other functions (this function will work with any metric).
#' @param n The number of desired ntiles. Default = 3.
#' @param labels Text labels for ntiles. Default = c("low", "medium", "high")
#' @param sampledata A data frame containing the sample annotation data. Must contain both sample names and the faceted group "PPBC".
#' @param survival_type Must be either "os" (overall survival) or "drs" (distant recurrence) 
#'
#' @return A list of two plots, one with the unadjusted survival curves by PPBC group, and one with the conditionally 
#' adjusted for clinical covariates survival curves for the entire dataset. See survminer::ggadjustedcurves for details.
#' @export
#'
#' @examples km_ntiles("MS4A1")
#' 
km_ntiles <- function(gene_id, id_type = "symbol", gene_dict = gx_annot, geneEx = ens_mat,
                      n = 3, labels = c("low", "medium", "high"),
                      sampledata = sample_data, survival_type = "os"){
  
  stopifnot(survival_type %in% c("os", "drs"))
  ntiles <- gene_ntile(gene_id = gene_id, id_type = id_type, gene_dict = gene_dict,
                       geneEx = geneEx, n = n, labels = labels)
  
  sd <- dplyr::left_join(ntiles, sampledata, by = "sample_name")
  #  return(sd)
  
  #Toggle for survival type
  if(survival_type == "os"){
    ti <- "months_overall_survival"
    ev <- "overall_survival"
    ylab <- "Overall survival probability"
    type <- "Overall survival"
  }
  
  if(survival_type == "drs"){
    ti <- "months_to_drs"
    ev <- "distant_recurrence"
    ylab <- "Distant recurrence probability"
    type <- "Distant recurrence"
  }
  
  covariates <- c("age_at_diagnosis", "year_of_diagnosis", "grade",
                  "stage", "surgery", "radiotherapy", "hormonetherapy",
                  "chemotherapy", "herceptin", "PAM50", "study_group")
  
  cov <- paste(covariates, collapse = "+")
  svf = paste0('Surv(time=',ti,', event=',ev,')~')
  facet_formula <- as.formula(paste0(svf, "labels"))
  reduced =  as.formula(paste0(svf, cov))
  formula = as.formula(paste0(svf, paste0(cov, "+ntile")))
  
  reslist <- list()
  gn <- unique(sd$gene_symbol)
  
  #Facet plot: Show unadjusted curves for gene by study group
  facet_plot <- survminer::ggsurvplot(
    #fit = survfit(Surv(months_overall_survival, overall_survival)~ labels, data = as.data.frame(sd)), 
    #If you want to use a formula object, use surv_fit instead of survfit
    fit = survminer::surv_fit(facet_formula, data = as.data.frame(sd)), 
    xlab = "Months", 
    ylab = ylab,
    title = paste(type, "with", gn, "expression by PPBC group"),
    pval = T, facet.by = "PPBC")
  
  reslist$facet_plot <- facet_plot
  
  #return(reslist)
  
  #Cox models for all covariates, study group and ntile gene expression
  c <- survival::coxph(formula, data = as.data.frame(sd))
  #print(c)
  
  #Store p value associated with the coefficient of the ntile
  fit_ntile <- as.data.frame(summary(c)$coefficients)
  p_ntile <- signif(fit_ntile[rownames(fit_ntile) == "ntile", ]$`Pr(>|z|)`, 2)
  
  #Anova p value vs reduced model
  d <- survival::coxph(reduced, data = as.data.frame(sd))
  aod <- anova(d,c) #Conventional to list from smallest to largest but order does not matter
  aod_p <- signif(aod[4][2,], 2)
  
  #Adjusted survival curve for all 4 PPBC groups
  #B too small to show groupwise
  #Use surv_adjustedcurves to get data and the then manually plot for more control over appearance
  adj_curv <- survminer::surv_adjustedcurves(fit = c, variable = "ntile", data = as.data.frame(sd), method = "conditional") %>%
    left_join(select(mutate(sd, ntile = as.factor(ntile)), ntile, labels), by = c("variable" = "ntile")) %>%
    ggplot(., aes(x = time, y = surv, color = labels)) + 
    geom_step(size = 1) + theme_survminer() + scale_y_continuous(limits = c(0,1)) + 
    ylab(label =  ylab) + labs(color = paste(gn, "expression")) +
    ggtitle(paste(type, "for gene", gn, "adjusted for clinical covariates")) +
    annotate(geom = "text", label = paste0("Coef pval for ntile: ", p_ntile), x= 25, y =0.1) +
    annotate(geom = "text", label = paste0("Anova pval vs reduced: ", aod_p), x= 25, y =0.2)
  
  reslist$adj_curv <- adj_curv
  
  return(reslist)
  
}

#' Plot Kaplan-Meier curves of ntiles of a given gene for one vs rest comparisons
#'
#' @param gene_id An ensembl id or gene symbol
#' @param id_type Must be either "symbol" or "ensembl"
#' @param gene_dict A data frame containing a column "gene name" with gene symbols and a column "ensembl_gene_id" with ensembl gene ids
#' @param geneEx A gene expression matrix with ensembl IDs as row names. Should be tmm/log2 normalized for compatibility with other functions 
#' (this function will work with any metric).
#' @param n The number of desired ntiles. Default = 3.
#' @param labels Text labels for ntiles. Default = c("low", "medium", "high")
#' @param sampledata A data frame containing the sample annotation data.
#'  Must contain both sample names and the column representing one vs rest comparison (1 and 0).
#' @param survival_type Must be either "os" (overall survival) or "drs" (distant recurrence)
#' @param return_list Whether to return a list of 4 plots instead of arranging them into one
#' @param p_method If "anova" (default), the p value is the full model including the gene ntile vs
#' a reduced model containing the clinical covariates only. If "coef", the p value is the
#' derived from the coeficient in the cox model corresponding to gene ntiles.
#'
#' @return Four plots arranged together: the unadjusted survival curves faceted by one vs rest group, the conditionally 
#' adjusted for clinical covariates survival curves for the group of interest, 
#' and the adjusted plot for the other groups. See survminer::ggadjustedcurves for details.
#' @export
#'
#' @examples km_ntiles("MS4A1")

km_ntiles_ovr <- function(gene_id, id_type = "symbol", gene_dict = gx_annot, geneEx = ens_mat,
                          n = 3, labels = c("low", "medium", "high"), ovr_column = "involution",
                          sampledata = sample_data, survival_type = "os",
                          p_method = "anova", return_list = F){
  
  stopifnot(survival_type %in% c("os", "drs"))
  stopifnot(p_method %in% c("anova", "coef"))
  stopifnot(ovr_column %in% colnames(sampledata))
  
  ntiles <- gene_ntile(gene_id = gene_id, id_type = id_type, gene_dict = gene_dict,
                       geneEx = geneEx, n = n, labels = labels)
  
  sd <- dplyr::left_join(ntiles, sampledata, by = "sample_name")
  #  return(sd)
  
  #Toggle for survival type
  if(survival_type == "os"){
    ti <- "months_overall_survival"
    ev <- "overall_survival"
    ylab <- "Overall survival probability"
    type <- "Overall survival"
  }
  
  if(survival_type == "drs"){
    ti <- "months_to_drs"
    ev <- "distant_recurrence"
    ylab <- "Distant recurrence probability"
    type <- "Distant recurrence"
  }
  
  covariates <- c("age_at_diagnosis", "year_of_diagnosis", "grade",
                  "stage", "surgery", "radiotherapy", "hormonetherapy",
                  "chemotherapy", "herceptin", "PAM50")
  
  cov_legend <- c("age", "year of diagnosis", "grade",
                  "stage", "treatment", "PAM50")
  
  cov <- paste(covariates, collapse = "+")
  svf = paste0('Surv(time=',ti,', event=',ev,')~')
  facet_formula <- as.formula(paste0(svf, "labels"))
  #ntile_formula <- as.formula(paste0(svf, "ntile"))
  reduced =  as.formula(paste0(svf, cov))
  formula = as.formula(paste0(svf, paste0(cov, "+ntile")))
  
  reslist <- list()
  gn <- unique(sd$gene_symbol)
  
  #Split sample data into those samples belonging to the group of interest (or not)
  sd1 = sd[sd[[ovr_column]] == 1,]
  sd0 = sd[sd[[ovr_column]] == 0,]
  
  
  curv1 <- survminer::ggsurvplot(
    fit = survminer::surv_fit(facet_formula, data = as.data.frame(sd1)), 
    xlab = "Months", 
    ylab = ylab,
    legend = "none",
    #title = paste(type, "(unadjusted), with", gn, "\nexpression in", ovr_column, "samples"),
    title = paste("Unadjusted curve for", ovr_column, "samples"),
    pval = T, #Logrank method
    )
  
  curv1 = curv1$plot #Otherwise we get a ggsurvplot object, which cannot be used with ggarrange
  
  curv1 <- curv1 +
    theme(plot.title = element_text(size = 12),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()
    )
  
  reslist$curv1 <- curv1
  
  #Unadajusted kaplan meier plot in other samples
  curv0 <- survminer::ggsurvplot(
    fit = survminer::surv_fit(facet_formula, data = as.data.frame(sd0)), 
    xlab = "Months", 
    ylab = ylab,
    legend = "none",
    #title = paste(type, "(unadjusted), with", gn, "\nexpression in other samples"),
    title = paste("Unadjusted curve for other samples"),
    pval = T, #Logrank method
    ) 
  
  curv0 <- curv0$plot #Otherwise we get a ggsurvplot object, which cannot be used with ggarrange
  
  curv0 <- curv0 +
    theme(plot.title = element_text(size = 12),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()
    )
  
  reslist$curv0 <- curv0
  
  #return(reslist)
  
  #Cox models for involution group, adjusted for all covariates and ntile gene expression
  c <- survival::coxph(formula, data = sd1)
  #print(c)
  
  #Store p value associated with the coefficient of the ntile
  fit_ntile <- as.data.frame(summary(c)$coefficients)
  p_ntile <- signif(fit_ntile[rownames(fit_ntile) == "ntile", ]$`Pr(>|z|)`, 2)
  if (p_ntile < 0.001){
    p_ntile <- "< 0.001"
  }
  
  #print(p_ntile)
  #Anova p value vs reduced model
  d <- survival::coxph(reduced, data = as.data.frame(sd1))
  aod <- anova(d,c) #Conventional to list models from smallest to largest but order does not matter
  aod_p <- signif(aod[4][2,], 2)
  if (aod_p < 0.001){
    aod_p <- "< 0.001"
  }
  
  if (p_method == "anova"){
    p <- aod_p
  } else {
    p <- p_ntile
  }
  
  #Adjusted survival curve for samples belonging to group of interest
  #Use surv_adjustedcurves to get data and the then manually plot for more control over appearance
  adj_curv1 <- survminer::surv_adjustedcurves(fit = c, variable = "ntile", data = as.data.frame(sd1), method = "conditional") %>%
    left_join(select(mutate(sd1, ntile = as.factor(ntile)), ntile, labels), by = c("variable" = "ntile")) %>%
    ggplot(., aes(x = time, y = surv, color = labels)) + 
    geom_step(size = 1) + theme_survminer() + scale_y_continuous(limits = c(0,1)) +
    #scale_x_continuous(limits = c(0,250)) +
    ylab(label =  ylab) +
    labs(color = paste(gn, "expression")) +
    #ggtitle(paste(type, "for gene", gn, "in", ovr_column, "\nsamples, adjusted for clinical covariates")) +
    ggtitle(paste("Adjusted curve for", ovr_column,"samples")) +
    #annotate(geom = "text", label = paste0("Coef pval for ntile: ", p_ntile), x= 25, y =0.1) +
    #annotate(geom = "text", label = paste0("Anova pval vs reduced: ", aod_p), x= 25, y =0.2) +
    annotate(geom = "text", label = paste("p =", p), x= 25, y =0.1) +
    theme(plot.title = element_text(size = 12),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "bottom"
    )
  
  reslist$adj_curv1 <- adj_curv1
  
  #Cox models for rest group, adjusted for all covariates and ntile gene expression
  c <- survival::coxph(formula, data = sd0)
  #print(c)
  
  #Store p value associated with the coefficient of the ntile
  fit_ntile <- as.data.frame(summary(c)$coefficients)
  p_ntile <- signif(fit_ntile[rownames(fit_ntile) == "ntile", ]$`Pr(>|z|)`, 2)
  if (p_ntile < 0.001){
    p_ntile <- "< 0.001"
  }
  
  #Anova p value vs reduced model
  d <- survival::coxph(reduced, data = as.data.frame(sd0))
  aod <- anova(d,c) #Conventional to list models from smallest to largest but order does not matter
  aod_p <- signif(aod[4][2,], 2)
  if (aod_p < 0.001){
    aod_p <- "< 0.001"
  }
  
  if (p_method == "anova"){
    p <- aod_p
  } else {
    p <- p_ntile
  }
  
  #Adjusted survival curve for samples belonging to group of interest
  #Use surv_adjustedcurves to get data and the then manually plot for more control over appearance
  adj_curv0 <- survminer::surv_adjustedcurves(fit = c, variable = "ntile", data = as.data.frame(sd0), method = "conditional") %>%
    left_join(select(mutate(sd0, ntile = as.factor(ntile)), ntile, labels), by = c("variable" = "ntile")) %>%
    ggplot(., aes(x = time, y = surv, color = labels)) + 
    geom_step(size = 1) + theme_survminer() + scale_y_continuous(limits = c(0,1)) + 
    #scale_x_continuous(limits = c(0,250)) +
    ylab(label =  ylab) +
    labs(color = paste(gn, "expression")) +
    #ggtitle(paste(type, "for gene", gn, "\nin other samples, adjusted for clinical covariates")) +
    ggtitle(paste("Adjusted curve for other samples")) +
    #annotate(geom = "text", label = paste0("Coef pval for ntile: ", p_ntile), x= 25, y =0.1) +
    #annotate(geom = "text", label = paste0("Anova pval vs reduced: ", aod_p), x= 25, y =0.2) +
    annotate(geom = "text", label = paste("p =", p), x= 25, y =0.1) +
    theme(plot.title = element_text(size = 12),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "bottom"
          )
  
  reslist$adj_curv0 <- adj_curv0
  
  if(return_list == T){
    return(reslist)  
  }
  
  # Arrange them all together
  arplot <- ggpubr::ggarrange(plotlist = reslist, nrow = 2, ncol = 2)
  
  arplot <- annotate_figure(arplot,
                  top = text_grob(paste(type, "for gene", gn), face = "bold", size = 14),
                  #bottom = text_grob("Time (Months)", size = 12),
                  bottom = text_grob(paste0("Time (Months)\n\nClinical covariates: ",
                                            paste(cov_legend, collapse=", ")),
                                     size = 12),
                  left = text_grob(ylab, rot = 90),
                  #right = "I'm done, thanks :-)!",
                  #fig.lab = "Figure 1", fig.lab.face = "bold"
  )
  
  return(arplot)
  
}

#' Combine and display results from survival analyses
#'
#' @param id An ensembl id or gene symbol
#' @param id_type Must be either "symbol" or "ensembl"
#' @param s A data frame containing overall survival results with the following columns: 
#' "gene_name,fdr,beta,HR (95% CI for HR),p.value,ensembl_gene_id,gene_type,description"
#' @param d A data frame containing distant recurrence results, with the same columns
#' @param ios A data frame containing overall survival results from the interaction model, containing at least the aforementioned columns.
#' May also contain other columns.
#' @param idrs A data frame containing distant recurrence results from the interaction model, containing at least the aforementioned columns.
#' May also contain other columns.
#' @param gene_dict A data frame containing a column "gene name" with gene symbols and a column "ensembl_gene_id" with ensembl gene ids

#' @return A data frame combining the above information
#' @export
#'
#' @examples
gene_survival <- function(id, id_type = "symbol", s = os, d = drs, ios = inv_int_os, idrs = inv_int_drs,
                          gene_dict = gx_annot){
  
  if(id_type == "ensembl"){
    ens_id <- id
  } else {
    gnt <- gene_dict[gene_dict$gene_name == id, "ensembl_gene_id"]
    ens_id <- unique(gnt$ensembl_gene_id)
    if (length(ens_id) > 1) {
      warning(paste0("Multiple ensembl ids retrieved for ", id, ": ",
                     paste(ens_id,collapse=","), ". Returning first in series."))
      ens_id <- ens_id[1]
    }
  }
  
  
  g_os <- s[s$ensembl_gene_id == ens_id, , drop=F]
  g_os <- g_os %>% mutate(type = "overall survival") %>%
    select(type, everything())
  
  g_drs <- d[d$ensembl_gene_id == ens_id, , drop=F]
  g_drs <- g_drs %>% mutate(type = "distant recurrence") %>%
    select(type, everything())
  
  g_int_os <- ios[ios$ensembl_gene_id == ens_id, , drop=F]
  g_int_os <- g_int_os %>% mutate(type = "interaction os") %>%
    select(type, gene_name, fdr, beta, `HR (95% CI for HR)`, p.value, ensembl_gene_id, gene_type, description)
  
  g_int_drs <- idrs[idrs$ensembl_gene_id == ens_id, , drop=F]
  g_int_drs <- g_int_drs %>% mutate(type = "interaction drs") %>%
    select(type, gene_name, fdr, beta, `HR (95% CI for HR)`, p.value, ensembl_gene_id, gene_type, description)
  
  g <- bind_rows(g_os, g_drs, g_int_os, g_int_drs)
  
  return(g)
}

#' Combine and display diffex results
#'
#' @param list_reports A list of results from DESeq2 in data frame format.
#' @param ensembl_id An ensembl gene id (character)
#' @param pthresh An FDR threshold below which to consider a gene significant
#' @param abslogfcthresh A log2fold change threshold below which a gene is significant
#'
#' @return A data frame summarizing the significance of a gene amongst multiple DESeq2 results
#' @export
#'
#' @examples
diffex_report = function(list_reports, ensembl_id, pthresh = 0.05, abslogfcthresh = 0.5){
  
  #print(paste("Adjusted p value threshold:", pthresh, ", abs(log2FoldChange) threshold:", abslogfcthresh))
  
  comp = lapply(list_reports,
                function(x) x[x$ensembl_gene_id==ensembl_id,
                              c("ensembl_gene_id","gene_name","log2FoldChange","padj")]) %>%
    unlist() %>% tibble::enframe() %>% tidyr::separate(name, into = c("comparison", "field"), sep="\\.") %>%
    tidyr::spread(key=field, value=value) %>%
    dplyr::mutate(log2FoldChange = as.numeric(log2FoldChange),
           padj = as.numeric(padj),
           padj_cutoff = !!pthresh,
           absLog2FoldChange_cutoff = !!abslogfcthresh)
  
  comp = comp %>% dplyr::mutate(sig = (padj < !!pthresh) & (abs(log2FoldChange) > !!abslogfcthresh)) %>%
    dplyr::select(comparison, sig, padj, log2FoldChange, everything())
  
  comp = comp %>% dplyr::arrange(padj)
  
  return(comp)
}

#' Combine survival and diffex results
#'
#' @param id An ensembl id or gene symbol
#' @param id_type Must be either "symbol" or "ensembl"
#' @param ... Other parameters to pass to gene_survival()
#'
#' @return A list of data frames containing both
#' @export
#'
#' @examples
gene_summary <- function(id, id_type = "symbol", ...){
  g <- gene_survival(id = id, id_type = id_type, ...)
  ens_id <- unique(g$ensembl_gene_id)
  r <- diffex_report(res_list, ens_id)
  res <- list(survival_report = g, diffex_report = r)
  return(res)
}


#' Make beehive and violin plots of genes from TMM/log2 transformed expression data
#'
#' @param id An ensembl id or gene symbol
#' @param id_type Must be either "symbol" or "ensembl"
#' @param ensembl_mat A matrix with TMM/log2 transformed expression data, ensembl ids as rownames and samples as columns
#' @param symbol_mat A matrix with TMM/log2 transformed expression data, gene names/symbols as rownames and samples as columns
#' @param sampledata A data frame containing the sample annotation data. Must contain both sample names and the faceted group "PPBC".
#'
#' @return A list of plots, beehive and violin
#' @export
#'
#' @examples
tmm_plots <- function(id, id_type = "symbol", ensembl_mat = ens_mat, symbol_mat = sym_mat, sampledata = sample_data){
  
  stopifnot(id_type %in% c("symbol", "ensembl"))
  
  if(id_type == "symbol"){
    mat <- symbol_mat
  } else {
    mat <- ensembl_mat
  }
  
  mat <- mat[rownames(mat) == id, , drop=F]
  
  if(nrow(mat) == 0){
    stop("Id not found in expression matrix, do you have the right identifier?")
  }
  
  df <- as.data.frame(mat) %>% rownames_to_column("gene_name") %>%
    gather(key = "sample_name", -gene_name, value = "tmm_log") %>%
    left_join(., sample_data, by = "sample_name")
  
  #return(df)
  
  bh <- df %>%
    ggplot(aes(x = factor(PPBC, levels = c("nulliparous", "pregnant", "lactation", "involution")), y = tmm_log)) +
    geom_jitter(aes(color = PAM50), height = 0, width = 0.2) +
    geom_boxplot(alpha = 0) +
    xlab("PPBC") +
    ylab("Log(TMM expression)") +
    ggthemes::scale_color_colorblind() +
    ggthemes::theme_clean() +
    ggtitle(paste("TMM/Log normalized expression of", unique(df$gene_name)))
  
  vb <- df %>%
    ggplot(aes(x = gene_name, y = tmm_log)) +
    geom_jitter(aes(color = study_group), height = 0, width = 0.2) +
    geom_violin(alpha = 0 ) +
    xlab("PPBC") +
    ylab("Log(TMM expression)") +
    ggthemes::scale_color_colorblind() +
    ggthemes::theme_clean()
  
  #rsl <- ggpubr::ggarrange(plotlist = list(bh, vb), ncol = 2)
  rsl <- list(beehive = bh, violin = vb)
  
  return(rsl)
}



#Not very useful, get rid of this?
cox_dotplot <- function(gene_id,  cox_results, id_type = "gene_name", cox_events = "os",
                        data = coxdata, gene_annotation = gx_annot){
  
  stopifnot(id_type %in% c("gene_name", "ensembl_gene_id"))
  stopifnot(cox_events %in% c("os", "drs")) 
  
  if(cox_events == "os"){
    event_data = coxdata %>% select(sample_name, study_group, months_overall_survival, overall_survival)
    title_string = "OS for interaction "
    legend_label = "Death"
  } else {
    event_data = coxdata %>% select(sample_name, study_group, months_to_drs, distant_recurrence)
    title_string = "DRS for interaction "
    legend_label = "Metastasis"
  }
  
  event_type = colnames(event_data)[3:4]
  colnames(event_data)[3:4] <- c("time_months", "event")
  
  if(id_type == "gene_name"){
    gene_annot <- gx_annot[gx_annot$gene_name==gene_id,,drop=F]
  } else {
    gene_annot <- gx_annot[gx_annot$ensembl_gene_id==gene_id,,drop=F]
  }
  
  #stopifnot(nrow(gene_annot) == 1)
  if (nrow(gene_annot) > 1){
    gene_annot <- head(gene_annot, 1)
  }
  
  gene_data = coxdata[,c(gene_annot$ensembl_gene_id), drop=F]
  gene_data = cbind(gene_data, event_data)
  
  colnames(gene_data)[1] = "TMM_lognorm_expression"
  gene_data$gene_name = gene_annot$gene_name
  
  gene_results = cox_results[cox_results$ensembl_gene_id == gene_annot$ensembl_gene_id, ]
  stopifnot(nrow(gene_results) == 1)
  
  title_string = paste0(title_string, paste0("(involution*", gene_annot$gene_name), ")\n") #Response type for gene name
  title_string = paste(title_string,
                       paste(paste(colnames(gene_results)[2:4], gene_results[2:4], sep=": "), #FDR, beta and HR
                             collapse = ", "))
  
  gene_data %>%
    mutate(event = as.factor(event),
           study_group = recode(study_group,
                                non_prbc = "Nulliparous",
                                ppbc_inv = "Involution",
                                ppbc_lac = "Lactation",
                                prbc = "Pregnancy",
                                .default = levels(study_group))) %>%
    ggplot(aes(x =time_months, y = TMM_lognorm_expression, color = event)) +
    geom_point() +
    geom_rug()+
    #geom_smooth(method='lm',se=F)+
    facet_wrap(~study_group)+
    scale_color_colorblind() + theme_few() +
    ggtitle(title_string) +
    ylab("log(TMM)") +
    xlab("Time (Months)") +
    labs(color=legend_label)
  
}