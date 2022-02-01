library(conflicted)

library(dplyr)
library(purrr)
library(spatstat)
library(yaml)

options(warn = 2)

segment_tissue_density <- function(cells, tissues, resolution, minimum_density) {
  # First build spatial point pattern with bounding box window such that we can use it with
  # spatstat.
  cells_pp <- with(cells, ppp(
    x=x,
    y=y,
    marks=cell_type,
    window=owin(xrange = c(min(x), max(x)), yrange = c(min(y), max(y)))))
  # Improve window with a convex boundary estimate.
  # FIXME: Use the actual boundaries
  Window(cells_pp) <- ripras(cells_pp)
  tissue_density <- map(tissues, function (a) {
    tissue_pp <- subset(cells_pp, marks %in% a$cell_types)
    density(tissue_pp, sigma=a$kernel_bandwidth, eps=resolution)
  })
  # Reduce list of tissue densities to the maximum density and the tissues with that maximum
  # density. Works with spatstat images to keep dimensions and units.
  max_density <- reduce2(tissue_density, names(tissue_density),
    function (accumulation, density, tissue) {
      higher_mask <- density > accumulation$density
      accumulation$density[higher_mask] <- density[higher_mask]
      accumulation$tissue[higher_mask] <- tissue
      accumulation
    },
    # Density needs to be higher than a minimum value to be assigned a tissue
    .init=list(
      density = as.im(minimum_density, tissue_density[[1]]),
      tissue = as.im(NA_character_, tissue_density[[1]])))
  max_density$tissue
}

if (sys.nframe() == 0) {
  parse_args <- function(args) {
    arg_names <- c('cells', 'config', 'out') 
    stopifnot(length(args) == length(arg_names))
    args <- as.list(args)
    names(args) <- arg_names
    config <- yaml.load_file(args$config)
    args <- append(args, config)
    args
  }

  args <- parse_args(commandArgs(T))

  cells <- readRDS(args$cells)
  segmentation <- segment_tissue_density(
    cells=cells,
    tissues=args$tissues,
    resolution=args$resolution,
    minimum_density=args$minimum_density)
  saveRDS(segmentation, args$out)
}
