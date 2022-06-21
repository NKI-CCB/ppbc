library(conflicted)

library(dplyr)
require(spatstat.core)

source('src/spatial/model_spatstat.R')

options(warn = 2)

pixel_size <- 0.0005 #mm
unitname <- c('milimeter', 'milimeters')


compute_Lcross <- function(cells, radii, cell_type1, cell_type2, ...) {
  if (sum(spatstat.geom::marks(cells) == cell_type2) == 0) {
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
  purrr::map_dfr(rlang::set_names(immune_cell_types), function (ct) {
      cells <- subset(cells, marks %in% c(ct, 'PanCK+'))
      with_simulation(cells, nsim, compute_Lcross, subset_poisson_null, radii = radii,
                      cell_type1 = 'PanCK+', cell_type2 = ct, resample_mark = ct,
                      .progress_bar = progress_bar)
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

  cells <- read_cells(args$objects, args$annotation)
  res <- model_Lcross(cells, args$progress_bars)
  readr::write_tsv(res, args$output)
}
