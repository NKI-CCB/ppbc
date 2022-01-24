library(here)
library(tidyverse)

countDir <- "results/spatial/marker_cell_counts"

#Summarize by sample and panel, include batch
sum_cell_counts <- function(inFiles = list.files(here(countDir))){
  lapply(inFiles, function(x){
    info <- tibble(filename = x)
    #handle "T20-62169_I1"  "T20-62169_II1"
    info <- info %>%
      mutate(filename = str_replace(filename, "_I", "-I"))
    info <- info %>%
      tidyr::separate(col = filename, into = c("t_number", "panel", "batch"),
                      sep = "_") %>%
      mutate(batch = str_remove(batch, ".csv"))
    
    df <- read_csv(here(file.path(countDir, x)), show_col_types = F)
  bind_cols(info, df)
}) %>%
  bind_rows() %>%
    mutate(classifier_label = str_replace(classifier_label, "tumor", "Tumor"))
}

x <- list.files(here(countDir))
fix <- x[str_detect(x, "_I|_II1")]


