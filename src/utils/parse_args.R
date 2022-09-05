library(dplyr)
library(stringr)

retrieve_args <- function(args = commandArgs(trailingOnly=TRUE)){
  #Is a vector with space as a delimiter
  print(args)
  
  stopifnot(length(args) != 0)
  #If every argument is named, will be divisible by two
  stopifnot(length(args) %% 2 == 0)
  
  #Names are
  argnames <- args[stringr::str_detect(args, "--")]
  stopifnot(stringr::str_sub(argnames, 1, 2) == '--')
  argvals <- args[!stringr::str_detect(args, "--")]
  
  argdf <- tibble(
    argname = argnames,
    argval = argvals
  ) %>%
    mutate(argname = stringr::str_remove(argname, "--"))
  
  argdf
}