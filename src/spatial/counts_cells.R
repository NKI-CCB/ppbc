library(conflicted)

library(dplyr)
requireNamespace('ncdf4')
library(purrr)
library(readr)
library(stringr)
library(tidyr)

source('src/spatial/read_cells.R')

options(warn = 2)

parse_args <- function(args) {
  arg_names <- c('objects', 'out_fn') 
  stopifnot(length(args) == length(arg_names))
  args <- as.list(args)
  names(args) <- arg_names
  args
}

read_and_count_objects <- function(fn) {
  ds <- ncdf4::nc_open(fn)
  cell_df <- nc_read_data_frame(ds, "cell")
  positive <- nc_read_matrix(ds, "positive_classification") > 0
  dye <- str_split_fixed(ncvar_get(ds, 'dye'), ' ', 2)[,1]
  nc_close(ds)
  cell_df$marker_pos = apply(positive, 2, function (positive_mask) {
    paste0(dye[positive_mask], '+', collapse='_')
  })
  cell_df %>%
    group_by(marker_pos, classifier_label, .drop=F) %>%
    summarize(n=n(), .groups="drop")
}

if (sys.nframe() == 0) {
  args <- parse_args(commandArgs(T))
  cell_counts <- read_and_count_objects(args$objects)
  write_csv(cell_counts, args$out_fn)
}
