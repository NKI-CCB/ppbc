library(conflicted)

library(dplyr)
library(stringr)

source("src/spatial/read_cells.R")

options(warn = 2)

disallowed_marker_pairs <- list(
  c("CD3", "CD20"),
  c("CD8", "CD20"),
  c("FoxP3", "CD20")
)


# dplyr::case_when variant that outputs factors with consistent levels
case_when_fct <- function(...) {
  levels <- purrr::map_chr(rlang::list2(...), rlang::f_name)
  factor(dplyr::case_when(...), levels=levels)
}

correct_double_positivity <- function(objects, marker1, marker2, verbose=F){
  if (verbose) {
    print(paste("Correcting for: ", paste0(marker1, "/", marker2), "double positivity"))
  }
  
  marker1_pos_col <- paste0(marker1, "_positive")
  marker2_pos_col <- paste0(marker2, "_positive")
  
  double_positive_idx <- which(objects[, marker1_pos_col] & objects[, marker2_pos_col])
  marker1_intensity <- objects[double_positive_idx, paste0(marker1, "_intensity")]
  marker2_intensity <- objects[double_positive_idx, paste0(marker2, "_intensity")]
 
  #Correct the double positive rows by marking the lower of the two markers as false
  marker1_higher <- marker1_intensity > marker2_intensity
  objects[double_positive_idx, marker1_pos_col] <- marker1_higher
  objects[double_positive_idx, marker2_pos_col] <- !marker1_higher

  objects
}

correct_marker_positivity <- function(objects) {
  for (dmp in disallowed_marker_pairs) {
    marker1_in_panel <- paste0(dmp[1], "_positive") %in% colnames(objects)
    marker2_in_panel <- paste0(dmp[2], "_positive") %in% colnames(objects)
    if (marker1_in_panel & marker2_in_panel) {
      objects <- correct_double_positivity(objects, dmp[1], dmp[2])
    }
  }
  objects
}

define_cell_types <- list(
  MPIF26 = function(objects) {
    mutate(objects, cell_type = case_when_fct(
      # If a cell is FoxP3 positive, we assume it must be CD3+
      FoxP3_positive ~ "FoxP3+",
      CD3_positive ~ "CD3+FoxP3-",
      CD20_positive ~ "CD20+",
      CD27_positive ~ "CD27+CD20-CD3-",
      PanCK_positive ~ "PanCK+",
      !CD27_positive & !FoxP3_positive & !CD3_positive & !CD20_positive
        & !PanCK_positive ~ "Other"),
      .after=object_id)
  },
  MPIF27 = function(objects) {
    mutate(objects, cell_type = case_when_fct(
      CD8_positive & CD3_positive ~ "CD3+CD8+",
      CD3_positive ~ "CD3+CD8-",
      CD8_positive ~ "CD3-CD8+",
      CD20_positive ~ "CD20+",
      PanCK_positive ~ "PanCK+",
      # CD138+ is ignored if it's not on a CD20+ cell
      !CD8_positive & !CD3_positive & !CD20_positive & !PanCK_positive ~ "Other"),
      .after=object_id)
  })

if (sys.nframe() == 0) {
  parse_args <- function(args) {
    arg_names <- c('in_fn', 'panel', 'out_fn')
    stopifnot(length(args) == length(arg_names))
    args <- as.list(args)
    names(args) <- arg_names
    args
  }

  args <- parse_args(commandArgs(T))

  objects <- read_cells(args$in_fn, intensity=T)

  # Assume that every maker is reported in a single location
  objects <- objects[,sapply(objects, function(x) all(!is.nan(x)))]
  colnames(objects) <- str_replace(colnames(objects), "_(nucleus|cytoplasm|membrane)_", "_")

  objects <- correct_marker_positivity(objects)
  objects <- define_cell_types[[args$panel]](objects)
  stopifnot(!is.na(objects$cell_type))
  saveRDS(objects, args$out_fn)
}
