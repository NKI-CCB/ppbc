library(conflicted)

library(dplyr)
library(ggplot2)
library(spatstat)

options(warn = 2)

plot_segmentation <- function(segmentation, title) {
    as.data.frame(segmentation) %>%
      ggplot(aes(x = x, y = y)) + 
      geom_raster(aes(fill=value)) +
      scale_fill_brewer("Segmentation class", palette = "Set1", drop = FALSE) +
      scale_x_continuous(paste0("(", unitname(segmentation), ")")) +
      scale_y_continuous(paste0("(", unitname(segmentation), ")")) +
      labs(title = title)
}

if (sys.nframe() == 0) {
  parse_args <- function(args) {
    stopifnot(length(args) >= 2)
    list(in_fns = args[1:(length(args)-1)], out_fn = args[length(args)])
  }

  args <- parse_args(commandArgs(T))
  pdf(args$out_fn)
  for (fn in args$in_fns) {
    segmentation <- readRDS(fn)
    plt <- plot_segmentation(segmentation, title=fs::path_file(fn))
    # Some images have gaps, making geom_raster emit a warning when printing. There is easy way
    # to suppress specific warnings, so we suppress all warnings.
    options(warn=0)
    print(plt)
    options(warn=2)
  }
  dev.off()
}

