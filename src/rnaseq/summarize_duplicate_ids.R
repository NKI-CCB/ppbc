summarize_expression_duplicate_ids <- function(mat, id_column, f=colMeans, final_gene_symbol_colname="GeneSymbol"){
  require(dplyr)
  
  print(paste("Starting with gene expression matrix containing", nrow(mat), "rows."))
  
  #Easiest way to write functions with dplyr is to standardize the column name
  
  if(id_column != "symbol"){
    colnames(mat)[colnames(mat)==id_column] <- "symbol"
  }
  
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