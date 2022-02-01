library(conflicted)

library(spatstat)

options(warn = 2)

if (sys.nframe() == 0) {
  parse_args <- function(args) {
    stopifnot(length(args) >= 2)
    list(in_fns = args[1:(length(args)-1)], out_fn = args[length(args)])
  }

  args <- parse_args(commandArgs(T))
  pdf(args$out_fn)
  for (fn in args$in_fns) {
    segmentation <- readRDS(fn)
    plot(segmentation, main=fs::path_file(fn))
  }
  dev.off()
}

