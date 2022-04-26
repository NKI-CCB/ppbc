library(conflicted)

library(dplyr)

options(warn = 2)

pixel_size <- 0.0005 #mm
unitname <- c('milimeter', 'milimeters')

read_cells <- function(objects_fn, annotation_fn) {
  objects <- readRDS(objects_fn)
  max_size <- 1e6
  annot_bin <- readBin(annotation_fn, 'raw', max_size)
  if (length(annot_bin) == max_size) {
    stop('Buffer not large enough') # Fixme: loop until everything has been read
  }
  ow <- spatstat.geom::as.owin(sf::st_as_sfc(annot_bin) * pixel_size)
  spatstat.geom::unitname(ow) <- unitname
  ow <- spatstat.geom::dilation(ow, pixel_size) # Otherwise not all cells lie exactly within
  samples_known_points_outside_window <- c(
    'T21-60592_MPIF27_batch3', 
    'T20-62178_MPIF27_batch7') # FIXME: This sample is missing an area in the annotation
  if (any(stringr::str_detect(objects_fn, samples_known_points_outside_window))) {
    objects <- objects[spatstat.geom::inside.owin(objects$x, objects$y, ow), ]
  }
  objects_pp <- spatstat.geom::ppp(x=objects$x, y=objects$y, marks=objects$cell_type, window = ow)
  objects_pp
}

with_poisson_null <- function(cells, nsim, fun, ..., .progress_bar=NULL) {
  progress_bar <- .progress_bar
  # Calculating observed function here to allow future expansion where the observed estimate
  # is compared to distribution under the null hypothesis, for example when calculating
  # significance.
  observed_res <- fun(cells, ...)
  if (!is.null(progress_bar)) progress_bar$tick()
  null_res <- purrr::map_dfr(rpois(nsim, spatstat.geom::npoints(cells)), function (n_cells_sim)  {
      cells_sim <- spatstat.core::rpoint(n_cells_sim, win=spatstat.geom::as.owin(cells))
      res <- fun(cells_sim, ...)
      if (!is.null(progress_bar)) progress_bar$tick()
      res
    }, .id='simulation') %>%
    group_by(across(-c(simulation, estimate))) %>%
    summarise(null_estimate = mean(estimate), .groups = 'drop')
  left_join(observed_res, null_res, by = setdiff(colnames(observed_res), 'estimate')) %>%
    mutate(
      observed_estimate = estimate,
      delta_estimate = estimate - null_estimate,
      model = 'Poisson')
}

compute_L <- function(cells, radii) {
  spatstat.core::Lest(cells, r=radii, correction='none') %>%
    as_tibble() %>%
    transmute(
      measure = 'L',
      radius = r,
      estimate = un)
}

compute_spatstats <- function(cells, progress_bars, nsim=100) {
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
      with_poisson_null(subset(cells, marks==ct), nsim, compute_L, radii = radii,
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

  cells <- read_cells(args$objects, args$annotation)
  res <- compute_spatstats(cells, args$progress_bars)
  readr::write_tsv(res, args$output)
}
