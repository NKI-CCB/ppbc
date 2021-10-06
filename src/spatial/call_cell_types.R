library(conflicted)

library(dplyr)
library(stringr)

options(warn = 2)

disallowed_marker_pairs <- list(
  c("CD3", "CD20"),
  c("CD8", "CD20"),
  c("FoxP3", "CD20")
)

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

if (sys.nframe() == 0) {
  parse_args <- function(args) {
    arg_names <- c('in_fn', 'out_fn') 
    stopifnot(length(args) == length(arg_names))
    args <- as.list(args)
    names(args) <- arg_names
    args
  }

  args <- parse_args(commandArgs(T))

  objects <- readRDS(args$in_fn)

  # Assume that every maker is reported in a single location
  objects <- objects[,sapply(objects, function(x) all(!is.nan(x)))]
  colnames(objects) <- str_replace(colnames(objects), "_(nucleus|cytoplasm|membrane)_", "_")

  objects <- correct_marker_positivity(objects)

  #TODO: include calling of cell types

  saveRDS(objects, args$out_fn)
}
