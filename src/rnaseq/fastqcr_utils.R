retrieve_fastqc_module = function(qc, sample_annot, module.type,
                                  dir=here("results/fastqc"),
                                  failorwarn = c("FAIL", "WARN"),
                                  min.lib.size = 10^6.7){
  
  qc <- qc %>% filter(module == !!module.type) %>%
    filter(status %in% failorwarn) %>%
    filter(as.integer(tot.seq) > min.lib.size)
  
  qc = left_join(qc, select(sample_annot, fastq,
                            patient_ID, sample_name),
                 by = c("sample"="fastq"))
  
  file.paths = file.path(dir, qc$fastq)
  
  or = fastqcr::qc_read_collection(file = file.paths,
                                   sample_names = qc$sample_name,
                                   modules = module.type,
                                   verbose=F)[[1]]
  or = left_join(or, sample_annot, by = c("sample"="sample_name"))
  
  or = left_join(or, select(qc, -module), by=c("sample_id"="sample"))
  
  return(or)
}

add_qc_to_metadata <- function(metadata, qcdata, module.type){
  
  qcdata <- qcdata %>%
    filter(module == !!module.type) %>%
    dplyr::rename(status.module=status) %>%
    transmute(fastq = paste0(sample, ".fastq.gz"))
  
  newdata <- left_join(metadata, select(qcdata, sample, status.module),
                       by="fastq")
  colnames(newdata)[length(colnames(newdata))] = str_to_lower(str_replace_all(module.type, " ", "_"))
  return(newdata)
}
