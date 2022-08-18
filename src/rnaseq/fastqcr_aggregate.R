library(here)
library(fastqcr)

qc <- fastqcr::qc_aggregate(here("results/rnaseq/fastqc"))

saveRDS(qc, here("data/rnaseq/interim/02_fastqcr.Rds"))
