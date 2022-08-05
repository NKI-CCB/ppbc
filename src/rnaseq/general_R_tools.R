#' Read excel tables
#'
#' @description Reads all excel tabs from a given file into a list
#'
#' @param path The path to the excel file
#'
#' @return A list of data frame containing the content of the excel tabs.
#' The list names are the same as the tab names.
#'
read_excel_tabs = function(path){
  tab_names = readxl::excel_sheets(path)
  ldf = lapply(tab_names, function(x) readxl::read_excel(path, sheet = x))
  names(ldf) = tab_names
  return(ldf)
}


#' Subset TxImport
#'
#' @description Subset a 
#'
#' @param txi A list of matrices from tximport
#' @param cols A vector of column names to be kept
#'
#' @return A list containing the same matrices as produced by tximport: abundance, counts, length
#' containing only the column names in `cols`
#'
#' @examples
subset_tximport <- function(txdb, cols){
  lapply(txi, function(x) if(is.matrix(x)) return(x[, cols]) else return(x))
}



#' Generate a unique ID from multiple replicates
#' 
#' @description Creates a unique identifier for a series of samples with the same identifier within a data frame
#' 
#' @param df A data frame
#' @param column_name The column within the data frame containing the identifiers to be made unique
#'
#' @return A data frame that contains a new column "id" that contains the original sample name and a replicate
#' number separated by an underscore.
#'
#' @examples
#' df <- data.frame(sample = c(rep("a", 3), rep("b", 2)))
#' unique_id_from_replicates(df, "sample")
unique_id_from_replicates = function(df, column_name){
  
  stopifnot(column_name %in% colnames(df))
  if(c("rep", "id") %in% colnames(df)){
    stop("The rep and id columns will be overwritten") #FIXME
  }
  
  #Programming with dplyr the easy way: Standardize the column name
  #colnames(df)[colnames(df)==column_name] = "id"
  df$id = df[,column_name]
  
  #Columns for to track duplication and a dummy replicate column
  df = df %>% dplyr::mutate(dup = duplicated(id), rep=0) 
  
  #Add one to the previous counter if the id is duplicated
  rep=c()
  for (i in 1:nrow(df)){
    if(df[i,"dup"]==T){
      rep = c(rep, rep[[i-1]] + 1)
    } else {
      rep = c(rep, 1)
    }
  }
  df$rep = rep
  
  #Concatenate the rep counter with the id and ditch the extraneous columns
  df = df %>% dplyr::mutate(id = paste(id, rep, sep="_"))
  
  df <- df %>% dplyr::select(-dup) %>% dplyr::select(id, everything())
  
  return(df)
}

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

normTMMlog2 <- function(object){
  object = calcNormFactors(object, method="TMM")
  object = cpm(object, log=T, normalized.lib.sizes=T)
  return(object)
}

#Returns a factor with levels in the same order as they were supplied to case_when
fct_case_when <- function(...) {
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- levels[!is.na(levels)]
  factor(dplyr::case_when(...), levels=levels)
}