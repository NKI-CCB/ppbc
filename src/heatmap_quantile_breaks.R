#From https://slowkow.com/notes/pheatmap-tutorial/#quantile-breaks
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#mat_breaks <- quantile_breaks(mat, n = 11)

#pheatmap(
#  mat               = mat,
#  color             = inferno(length(mat_breaks) - 1),
# OR
#  color              = colorRampPalette(rev(brewer.pal(n = 11,
                                        #name = "RdYlBu")))(length(mat_breaks) - 1),
#  breaks            = mat_breaks,
#  border_color      = NA,
#  show_colnames     = FALSE,
#  show_rownames     = FALSE,
#  annotation_col    = mat_col,
#  annotation_colors = mat_colors,
#  drop_levels       = TRUE,
#  fontsize          = 14,
#  main              = "Quantile Color Scale"
#)

#Example heatmap with quantile normalization
library(here)
library(tidyverse)
library(pheatmap)

coxdata = readRDS("data/Rds/12_multi_genewise_os.Rd")

ann_df = coxdata[, c("sample_name","PPBC", "overall_survival")]
colnames(ann_df)[3] = "death"
ann_df = mutate(ann_df, death = as.factor(death))
ann_df = column_to_rownames(ann_df, "sample_name")

#Sort by survival
ann_df = ann_df[order(ann_df$death),]

#Subset geneEx to genes below threshold
if(nrow(filter(res, fdr <= !!pthresh)) < 5){
  mat = geneEx[rownames(geneEx) %in% head(res, 10)$gene_name, , drop=F]
  ttext = paste0("Top 10:\n ", analysis_type)
} else {
  mat = geneEx[rownames(geneEx) %in% filter(res, fdr <= !!pthresh)$gene_name, , drop=F] 
  ttext = paste0("Genewise fdr <", pthresh, ":\n ", analysis_type)
}
#dim(mat)

#Match matrix columns to annotation
mat = mat[,colnames(mat)[match(rownames(ann_df), colnames(mat))]]

#sp$overall_survival #white is too light
death_cols = sp$overall_survival
death_cols[2] = "lightgray"

#Quantile normalization
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(mat, n = 11)



pheatmap::pheatmap(mat, show_colnames = F,
                   scale="row", annotation_col = ann_df,
                   cluster_cols = T, drop_levels = T,
                   breaks =  mat_breaks,
                   color = colorRampPalette(rev(brewer.pal(n = 11,
                                                           name = "RdYlBu")))(length(mat_breaks) - 1),
                   annotation_colors = list(PPBC = ppbc_colors[-5],
                                            death = death_cols), #sp$overall_survival white is too light
                   main = ttext
) 