retrieve_fastqc_module = function(qc, sample_annot, module.type,
                                  dir=here("results/rnaseq/fastqc"),
                                  failorwarn = c("FAIL", "WARN"),
                                  min.lib.size = min_lib_size){
  
  qc <- qc %>% filter(module == !!module.type) %>%
    filter(status %in% failorwarn) %>%
    filter(as.integer(tot.seq) > min.lib.size)
  
  file.paths <- file.path(dir, qc$sample)
  
  or <- fastqcr::qc_read_collection(file = file.paths,
                                    sample_names = qc$sample,
                                    modules = module.type,
                                    verbose=F)$sequence_duplication_levels
  or <- left_join(mutate(or, fastq = paste0(sample, ".fastq.gz")),
                  sample_annot, by = c("fastq"))
  
  # or <- left_join(or, select(qc, -module), by=c("sample_id"="sample"))
  
  return(or)
}

add_qc_to_metadata <- function(metadata, qcdata, module.type){
  
  # Select only the status pertaining to the relevant module
  qcdata <- qcdata %>%
    filter(module == !!module.type) %>%
    dplyr::rename(status.module=status) %>%
    mutate(sample = paste0(sample, ".fastq.gz")) %>%
    rename(fastq=sample)
  
  # Add the module status to the metadata
  newdata <- left_join(metadata, select(qcdata, fastq, status.module),
                       by="fastq")
  
  #Rename the generic "status.module" column after the module it refers to
  colnames(newdata)[colnames(newdata)=="status.module"] <- str_to_lower(
    str_replace_all(module.type, " ", "_")) #prefer underscores over spaces
  
  newdata
}
