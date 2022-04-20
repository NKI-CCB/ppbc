library(conflicted)

library(dplyr)
library(purrr)
library(spatstat)
library(yaml)

options(warn = 2)

as_pp_cells <- function (cells, unitname = c("milimeter", "milimeters")) {
  # First build spatial point pattern with bounding box window such that we can use it with
  # spatstat.
  # FIXME: Use the actual boundaries from HALO
  cells_bb <- with(cells, owin(xrange = c(min(x), max(x)), yrange = c(min(y), max(y)),
                   unitname = unitname))
  cells_pp <- with(cells, ppp(
    x=x,
    y=y,
    marks=cell_type,
    unitname = unitname,
    window=cells_bb))
  # Improve window with a convex boundary estimate.
  Window(cells_pp) <- ripras(cells_pp)
  cells_pp
}

estimate_densities <- function(cells_pp, tissues, resolution) {
  map(tissues, function (a) {
    tissue_pp <- subset(cells_pp, marks %in% a$cell_types)
    density(tissue_pp, sigma=a$kernel_bandwidth, eps=resolution)
  })
}

segment_densities <- function(tissue_density, tissues, minimum_density) {
  tissue_wdensity <- map2(tissue_density, tissues, function(d, t) {d * t$weight})
  # Reduce list of tissue densities to the maximum density and the tissues with that maximum
  # density. Works with spatstat images to keep dimensions and units.
  max_density <- reduce2(tissue_wdensity, names(tissue_wdensity),
    function (accumulation, density, tissue) {
      higher_mask <- density > accumulation$density
      accumulation$density[higher_mask] <- density[higher_mask]
      accumulation$tissue[higher_mask] <- tissue
      accumulation
    },
    # Density needs to be higher than a minimum value to be assigned a tissue
    .init=list(
      density = as.im(minimum_density, tissue_wdensity[[1]]),
      tissue = as.im(factor(NA_character_, levels=names(tissue_wdensity)), tissue_wdensity[[1]])))
  max_density$tissue
}

segment_tissue_density <- function(cells, tissues, resolution, minimum_density) {
  cells_pp <- as_pp_cells(cells)
  tissue_density <- estimate_densities(cells_pp, tissues, resolution)
  segment_densities(tissue_density, tissues, minimum_density)
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
