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
    mutate(
      t_number = factor(t_number),
      batch = factor(batch),
    ) %>%
    group_by(classifier_label, t_number, panel, cell_type, batch, .drop=F) %>%
    summarize(
      area = ensure_one_value(classifier_area),
      n = n_distinct(cell)) %>%
    group_by(classifier_label) %>%
    mutate(
      area = ensure_one_value(area), # Fill in area when n==0
    ) %>%
    mutate(density = if_else(n==0, 0, n / area)) # Edge case if there are no cells at all
}

#Some aggregate groups for density (currently just for CD27+/- B and T cells)
agg_cells <- list(
  "CD3+FoxP3-" = c("CD3+CD27+FoxP3-", "CD3+CD27-FoxP3-"),
  "CD20+" = c("CD20+CD27+", "CD20+CD27-")
)

#Calculate density for a single aggregate group
density_aggs <- function(density, aggs){
  newgroup <- names(aggs)
  stopifnot(length(newgroup) == 1)
  oldgroups <- aggs[[1]]
  stopifnot(length(oldgroups) > 1)
  
  # Handle an edge case in which a sample is present is only one panel
  # (and not the one that contains the cell groups to be aggregated)
  if(!all(sapply(oldgroups, function(x){x %in% unique(density$cell_type)}))){
    return(density)
  }
  
  df <- density %>%
     dplyr::filter(cell_type %in% {{oldgroups}})
  
  df %>%
    group_by(t_number, classifier_label, panel, batch, area) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    mutate(density = if_else(n==0, 0, n / area), .after = n) %>%
    mutate(cell_type = {{newgroup}}, .after = panel)
}

#Loop over all aggregate groups
aggregate_densities <- function(density) {
  #Avoid duplicating entries if one row is absent
  res <- 0
  lapply(names(agg_cells), function(x){
    if(!all(sapply(agg_cells[x], function(y){y %in% unique(density$cell_type)}))){
      paste(paste0(unlist(agg_cells[[x]]), collapse = " and "), "not found")
    } else {
      res <- lapply(names(agg_cells), function(groups){
        density_aggs(density, aggs = agg_cells[groups])
      }) %>%
        #return aggregates
        bind_rows() %>%
        #recombine with input
        bind_rows(density, .)  
    }
  })
  
  if(res == 0){res <- density}
  res
}

if (sys.nframe() == 0) {
  args <- parse_args(commandArgs(T))
  objects <- readRDS(args$objects)
  density <- model_density(objects)
  all_densities <- aggregate_densities(density)
  write_tsv(all_densities, args$out)
}
