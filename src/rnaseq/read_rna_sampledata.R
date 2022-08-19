library(tidyverse)

#' Read RNA sampledata
#' 
#' @description Parse RNA sampledata. Duplicate patient IDs may exist where
#' a sample has been sequenced more than once. Column typing is therefore not 
#' set until later.
#'
#' @param path Character, the path to the tsv sampledata file
#' @param verbose Logical, whether to print column changes
#'
#' @return A tibble
#' @export
#'
#' @examples  read_rna_sampledata(sampledata_path)
read_rna_sampledata <- function(path, verbose = T){
  
  sampledata <- read_tsv(path, show_col_types = F) %>%
    filter(sample_type == "RNA", !is.na(sample_ID)) %>%
    rename(fastq=sample_ID)
  stopifnot(all(str_detect(sampledata$fastq, "fastq.gz")))
  
  if(verbose){print("Renaming sample_ID to fastq")}
  
  colnames_before <- colnames(sampledata)
  
  # Only include comments if there are any
  if(all(is.na(sampledata$comments))){
    sampledata <- sampledata %>% select(-comments)
  }
  
  # Experimental platform is redundant in this case
  # Batch HALO is not applicable
  # Ignore columns we don't know about
  stopifnot(all(sampledata$sample_type == "RNA"))
  stopifnot(all(sampledata$experimental_platform == "RNAseq"))
  if(verbose){print("Renaming batch_Leuven to batch")}
  sampledata <- sampledata %>%
    select(patient_ID, fastq, batch = batch_Leuven)
  
  if(verbose){
    print(paste("Columns kept:", paste(colnames(sampledata), collapse = ", ")))
    print(paste("Columns discarded:",
                paste(setdiff(colnames_before, colnames(sampledata)),
                      collapse = ", ")))
  }
  
  sampledata
}