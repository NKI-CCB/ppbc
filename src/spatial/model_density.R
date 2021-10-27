library(conflicted)

library(dplyr)
library(purrr)
library(readr)
library(tidyr)

options(
  warn = 2,
  dplyr.summarise.inform=F)

parse_args <- function(args) {
  arg_names <- c('objects', 'out') 
  stopifnot(length(args) == length(arg_names))
  args <- as.list(args)
  names(args) <- arg_names
  args
}

ensure_one_value <- function(x) {
  x <- unique(x)
  stopifnot(length(x) == 1)
  x
}

model_density <- function(objects) {
  objects %>%
    group_by(classifier_label, t_number, panel, cell_type) %>%  # FIXME: deal with batches
    summarize(
      area = ensure_one_value(classifier_area),
      n = n_distinct(cell)) %>%
    mutate(density = n / area)
}

if (sys.nframe() == 0) {
  args <- parse_args(commandArgs(T))
  objects <- readRDS(args$objects)
  density <- model_density(objects)
  write_tsv(density, args$out)
}
