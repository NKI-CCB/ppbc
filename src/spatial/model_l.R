library(conflicted)

library(dplyr)
require(spatstat.core)

source('src/spatial/model_spatstat.R')

options(warn = 2)

pixel_size <- 0.0005 #mm
unitname <- c('milimeter', 'milimeters')

#' Wrapper around spatstat.core::Lest that returns a tibble.
#'
#' @description Computes `spatstat.core::Lcross` on `cells` for `radii`. No
#'   correction is performed, this function is intended to be used in the `with_simulation` function
#'   from ``src/spatial/model_spatstat.R`  which should correct for the observation window.
#'
#' @param cells       A spatstat::ppp object with cells in an observation window.
#' @param radii       Radii to calculate L at, passed to the Lest function as the `r` argument.
#' @return A tibble with a `measure` column with 'L', a `r` column with the radius,
#' and an `estimate` column with the L  estimate.
compute_L <- function(cells, radii, ...) {
  spatstat.core::Lest(cells, r=radii, correction='none') %>%
    as_tibble() %>%
    transmute(
      measure = 'L',
      radius = r,
      estimate = un)
}

#' Calculate L of immune cells per immune cell type.
#'
#' @description Calculate the L based on distances between immune cells split by immune
#'   cell type. The marks on `cells` should be cell types, all marks that are not 'PanCK+' or
#'   'Other' are assumed to be immune cell types.
#'
#' @param cells       A spatstat::ppp object with cells in an observation window marked with
#'   cell types.
#' @param progress_bars A boolean to print progress bars to the terminal or not.
#' @param nsim        The number of simulations to run.
#' @return A tibble with results from `compute_L` adjusted for a null distribution that randomizes
#'   all cells. The tibble contains rows for different immune cell types and radii.
model_L <- function(cells, progress_bars, nsim=100) {
  cell_types <- levels(cells$marks)
  immune_cell_types <- setdiff(cell_types, c('PanCK+', 'Other'))

  if (progress_bars) {
    n_tick_progress = (nsim + 1) * length(immune_cell_types)
    progress_bar <- progress::progress_bar$new(
      format = "  simulating [:bar] :percent in :elapsed", total = n_tick_progress, clear = FALSE,
      width = 80)
  } else {
    progress_bar <- NULL
  }

  # Maximum radius to compute this statistic for is limited by the largest distance in our image
  # and what we are interested in.
  radius_max <- min(spatstat.geom::inradius(spatstat.geom::as.owin(cells)), 1.0)
  # Set the radii to compute, such they can later be aligned over multiple samples.
  radius_step <- 0.001
  radii <- seq(0, radius_max, radius_step)
  purrr::map_dfr(rlang::set_names(immune_cell_types), function (ct) {
      cells <- subset(cells, marks == ct)
      with_simulation(cells, nsim, compute_L, simple_poisson_null, radii = radii,
                      .progress_bar = progress_bar)
    }, .id = 'cell_type')
}

if (sys.nframe() == 0) {
  parse_args <- function(args) {
    arg_names <- c('objects', 'annotation', 'output', 'progress_bars') 
    stopifnot(length(args) == length(arg_names))
    args <- as.list(args)
    names(args) <- arg_names
    args$progress_bars <- as.logical(args$progress_bars)
    args
  }

  args <- parse_args(commandArgs(T))
  set.seed(123)
  cells <- read_cells(args$objects, args$annotation)
  res <- model_L(cells, args$progress_bars)
  readr::write_tsv(res, args$output)
}
