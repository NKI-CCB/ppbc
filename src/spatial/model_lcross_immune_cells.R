library(dplyr)
require(spatstat.core)

source('src/spatial/model_spatstat.R')

options(warn = 2)

pixel_size <- 0.0005 #mm
unitname <- c('milimeter', 'milimeters')


#' Wrapper around spatstat.core::Lcross that returns a tibble.
#'
#' @description Computes `spatstat.core::Lcross` on `cell_type1` and `cell_type2` for `radii`. No
#'   correction is performed, this function is intended to be used in the `with_simulation` function
#'   from ``src/spatial/model_spatstat.R`  which should correct for the observation window.
#'
#' @param cells       A spatstat::ppp object with cells in an observation window marked with
#'   cell types.
#' @param radii       Radii to calculate Lcross at, passed to that function as the `r` argument.
#' @param cell_type1   The cell type from which Lcross is measured, passed to Lcross as the second
#'   argument.
#' @param cell_type2   The cell type to which Lcross is measured, passed to Lcross as the third
#'   argument.
#' @return A tibble with cell types in `cell_type1` and `cell_type2`, a `measure` column with
#'   'Lcross', a `r` column with the radius, and an `estimate` column with the Lcross estimate.
#'   `estimate` columns is NA if no cells of cell_type1 or cell_type2 exist.
compute_Lcross <- function(cells, radii, cell_type1, cell_type2, ...) {
  # If the immune cell type has no cells (possible in some simulations), return NA.
  if (sum(spatstat.geom::marks(cells) == cell_type1) == 0
      | sum(spatstat.geom::marks(cells) == cell_type2) == 0) {
    tibble(estimate = NA)
  } else {
    spatstat.core::Lcross(cells, cell_type1, cell_type2, r=radii, correction='none') %>%
    as_tibble()
  } %>%
    transmute(
      cell_type1 = cell_type1,
      cell_type2 = cell_type2,
      measure = 'Lcross',
      radius = r,
      estimate = un)
}

#' Calculate L-cross between immune cells of different immune cell types.
#'
#' @description Calculate the L-cross based on distances from immune cells of one type to immune
#' cells of another type. The marks on `cells` should be cell types, all marks that are not 'PanCK+'
#' or 'Other' are assumed to be immune cell types.
#'
#' @param cells       A spatstat::ppp object with cells in an observation window marked with
#'   cell types.
#' @param progress_bars A boolean to print progress bars to the terminal or not.
#' @param nsim        The number of simulations to run.
#' @return A tibble with results from `compute_Lcross` adjusted a null distribution that randomized
#'   the immune cells of the first type but keeps the other cells fixed. The tibble contains rows
#'   for different immune cell types and radii.
model_Lcross <- function(cells, progress_bars, nsim=100) {
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
  radius_max <- min(spatstat.geom::inradius(spatstat.geom::as.owin(cells)), 0.1)
  # Set the radii to compute, such they can later be aligned over multiple samples.
  radius_step <- 0.001
  radii <- seq(0, radius_max, radius_step)
  purrr::cross2(cell_types, cell_types, .filter = `==`) %>%
    purrr::map_dfr(function (ctl) {
      cells <- subset(cells, marks %in% ctl)
      with_simulation(cells, nsim, compute_Lcross, subset_poisson_null, radii = radii,
                      cell_type1 = ctl[[1]], cell_type2 = ctl[[2]],
                      resample_mark = ctl[[2]], .progress_bar = progress_bar)
      })
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
  res <- model_Lcross(cells, args$progress_bars)
  readr::write_tsv(res, args$output)
}
