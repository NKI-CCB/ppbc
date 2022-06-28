# Shared functions for the different spatstat models

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

with_simulation <- function(cells, nsim, fun, simulation, ..., .progress_bar=NULL) {
  progress_bar <- .progress_bar
  # Calculating observed function here to allow future expansion where the observed estimate
  # is compared to distribution under the null hypothesis, for example when calculating
  # significance.
  observed_res <- fun(cells, ...)
  if (!is.null(progress_bar)) progress_bar$tick()
  null_res <- purrr::map_dfr(1:nsim, function (isim)  {
      cells_sim <- simulation$fun(cells, ...)
      res <- fun(cells_sim, ...)
      if (!is.null(progress_bar)) progress_bar$tick()
      res
    }, .id='simulation') %>%
    group_by(across(-c(simulation, estimate))) %>%
    summarise(null_estimate = mean(estimate, na.rm=T), .groups = 'drop')
  left_join(observed_res, null_res, by = setdiff(colnames(observed_res), 'estimate')) %>%
    mutate(
      observed_estimate = estimate,
      delta_estimate = estimate - null_estimate,
      model = simulation$name(...)) %>%
    select(-estimate)
}

simple_poisson_null <- list(
  fun = function(cells, ...) {
    n_cells <- spatstat.geom::npoints(cells)
    ow <- spatstat.geom::as.owin(cells)
    withCallingHandlers(
      spatstat.core::rpoint(n_cells, win=ow),
      # Ignore warning of creating a large number of points
      warning = function (cond) {
        if (stringr::str_starts(cond$message, "Attempting to generate ")) {
          invokeRestart("muffleWarning")
        }
    })
  },
  name = function (...) 'Poisson')

subset_poisson_null <- list(
  fun = function(cells, resample_mark, ...) {
    sample_cells <- subset(cells, marks == resample_mark)
    other_cells <- subset(cells, marks != resample_mark)
    n_cells <- spatstat.geom::npoints(sample_cells)
    ow <- spatstat.geom::as.owin(cells)
    simulated_cells <- withCallingHandlers(
      spatstat.core::rpoint(n_cells, win=ow),
      # Ignore warning of creating a large number of points
      warning = function (cond) {
        if (stringr::str_starts(cond$message, "Attempting to generate ")) {
          invokeRestart("muffleWarning")
        }
    })
    spatstat.geom::marks(simulated_cells) <- factor(
      resample_mark,
      levels = levels(spatstat.geom::marks(cells)))
    spatstat.geom::superimpose(simulated_cells, other_cells)
  },
  name = function(resample_mark, ...) {paste0("Poisson (", resample_mark, ")")})
