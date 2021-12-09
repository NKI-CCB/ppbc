library(conflicted)

library(dplyr)
library(readr)
library(purrr)

options(warn = 2)

if (sys.nframe() == 0) {
  # Read all tsv files on the command line and concatenate them into one big tsv
  parse_args <- function(args) {
    args <- commandArgs(T)
    stopifnot(length(args) >= 3)
    list(
      in_fns = args[1:(length(args) - 1)],
      out_fn = args[length(args)]
    )
  }

  args <- parse_args(commandArgs(T))

  df <- rlang::set_names(args$in_fns) %>%
    purrr::map_dfr(read_tsv)

  write_tsv(df, args$out_fn)
}
