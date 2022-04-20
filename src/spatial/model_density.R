library(conflicted)

library(dplyr)
library(purrr)
library(readr)
library(stringr)
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
  
  #Fix incidental lowercasing
  objects <- objects %>%
    mutate(classifier_label = stringr::str_to_sentence(classifier_label))
  
  region <- objects %>%
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
  
  # Also perform total aggregation on combined tumor and stroma
  totals <- region %>%
    group_by(t_number, panel, cell_type, batch) %>%
    summarise(area = sum(area), n = sum(n), .groups = "drop") %>%
    mutate(classifier_label = "Total", .before = everything()) %>%
    group_by(classifier_label) %>%
    mutate(
      area = ensure_one_value(area), # Fill in area when n==0
      classifier_label = factor(classifier_label, levels = c("Stroma", "Tumor", "Total"))
    ) %>%
    mutate(density = if_else(n==0, 0, n / area))
  
  
  bind_rows(region, totals)
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
    stop("This cell type is not found in this panel")
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
  
  #check to see whether the cells to be aggregated are present
  cell_presence <- sapply(agg_cells, function(x){x %in% unique(density$cell_type)})
  groups_kept <- colnames(cell_presence)[apply(cell_presence, 1, function(x){sum(x) == length(x)})]
  used_agg_cells <- agg_cells[names(agg_cells) %in% groups_kept]
  
  if(length(used_agg_cells)==0){return(density)}
  
  #Avoid duplicating entries if one row is absent
  res <- lapply(names(used_agg_cells), function(x){
      density_aggs(density, aggs = used_agg_cells[x])
    }) %>%
    #return aggregates
    bind_rows() %>%
  #recombine with input
  bind_rows(density, .) 
  
  #if(res == 0){res <- density}
  res %>% arrange(classifier_label, cell_type)
}

if (sys.nframe() == 0) {
  args <- parse_args(commandArgs(T))
  objects <- readRDS(args$objects)
  density <- model_density(objects)
  all_densities <- aggregate_densities(density)
  write_tsv(all_densities, args$out)
}


