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



