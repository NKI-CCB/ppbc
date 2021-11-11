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
  x <- x[!is.na(x)] # Because we store data in a tidy format, some values can disappear
  x <- unique(x)
  if (length(x) == 0) {
    NA
  } else if (length(x) == 1) {
    x
  } else {
    stop('Number of unique non-NA values > 1')
  }
}

model_density <- function(objects) {
  objects %>%
    group_by(classifier_label, t_number, panel, cell_type, .drop=F) %>%  # FIXME: deal with batches
    summarize(
      area = ensure_one_value(classifier_area),
      n = n_distinct(cell)) %>%
    mutate(density = if_else(n==0, 0, n / area)) # If n==0, we don't have area
}

if (sys.nframe() == 0) {
  args <- parse_args(commandArgs(T))
  objects <- readRDS(args$objects)
  density <- model_density(objects)
  write_tsv(density, args$out)
}
