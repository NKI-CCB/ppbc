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