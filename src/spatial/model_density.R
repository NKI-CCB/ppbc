library(conflicted)

library(dplyr)
library(purrr)
library(readr)
library(tidyr)

source('src/spatial/read_cells.R')

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

model_density <- function(fn) {
    cells_ds <- ncdf4::nc_open(fn)
    on.exit(ncdf4::nc_close(cells_ds), add=T)
    cells <- read_cells(cells_ds)
    cells <- left_join(
      cells,
      cells %>%
        pivot_longer(
          ends_with('positive'),
          names_to='dye',
          names_pattern='(.*)_positive',
          values_to='positive') %>%
        dplyr::filter(positive) %>%
        group_by(cell) %>%
        summarize(
          cell_type=factor(paste0(dye, '+', collapse='_'))),
      by='cell')
    area <- nc_read_matrix(cells_ds, 'area') %>%
      as_tibble() %>%
      pivot_longer(everything(), names_to='classifier_label', values_to='area')
    counts <- cells %>% group_by(cell_type, classifier_label) %>% summarize(n=n())
    left_join(counts, area, by='classifier_label')
}

if (sys.nframe() == 0) {
  args <- parse_args(commandArgs(T))
  density <- model_density(args$objects)
  density$density = density$n / density$area
  write_tsv(density, args$out)
}
