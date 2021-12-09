library('dplyr')
library(ncdf4)
requireNamespace('purrr')
requireNamespace('stringr')
library(tibble)
requireNamespace('tidyr')

pixel_size <- 0.0005  # mm
unitname <- c('milimeter', 'milimeters')

# Read a flag variable from a NetCDF dataset
nc_flagvar_get <- function(ds, var_name) {
    flag_values <- ncatt_get(ds, var_name, "flag_values")
    stopifnot(flag_values$hasatt)
    flag_meanings <- ncatt_get(ds, var_name, "flag_meanings")
    stopifnot(flag_meanings$hasatt)
    vals <- ncvar_get(ds, var_name)
    factor(vals, levels = flag_values$value, labels = flag_meanings$value)
}

# Read a data frame from a NetCDF4 dataset
nc_read_data_frame <- function(ds, dim) {
    v <- ds$dim[[dim]]$vals
    df <- tibble(v)
    names(df) <- c(dim)

    for (var_name in names(ds[["var"]])) {
        var <- ds[["var"]][[var_name]]
        if ((length(var[["dim"]]) == 1) &&
                (var[["dim"]][[1]][["name"]] == dim)) {
            if (ncatt_get(ds, var_name, "flag_meanings")$hasatt) { # Flag variable
                vals <- nc_flagvar_get(ds, var_name)
            } else {
                vals <- ncvar_get(ds, var_name)
            }
            df[[var_name]] <- vals
        }
    }
    df
}

# Read a matrix from a NetCDF4 dataset
nc_read_matrix <- function(ds, var_name) {
    mat <- ncvar_get(ds, var_name, collapse_degen=FALSE)
    row_dim <- ds[["var"]][[var_name]][["dim"]][[1]]
    col_dim <- ds[["var"]][[var_name]][["dim"]][[2]]
    mat_dim <- list(c(row_dim[["vals"]]), c(col_dim[["vals"]]))
    names(mat_dim) <- c(row_dim[["name"]], col_dim[["name"]])
    dimnames(mat) <- mat_dim
    mat
}

read_cells <- function(fn, intensity=F) {
    if (is.character(fn)) {
      ds <- ncdf4::nc_open(fn)
    } else {
      ds <- fn
    }
    cell_df <- nc_read_data_frame(ds, "cell")
    cell_matrices <- list()
    cell_matrices[['positive']] <- t(nc_read_matrix(ds, "positive_classification") > 0)
    if (intensity) {
      for (m in c('nucleus_intensity', 'cytoplasm_intensity', 'membrane_intensity')) {
        cell_matrices[[m]] <- t(nc_read_matrix(ds, m))
        stopifnot(nrow(cell_matrices[[m]]) == nrow(cell_df))
      }

    }
    cell_df <- cell_df %>%
        dplyr::mutate(
            x = pixel_size * (xmin + xmax) / 2,
            y = pixel_size * (ymin + ymax) / 2) %>%
        dplyr::distinct() %>%
        dplyr::rename_with(~ stringr::str_replace(., stringr::fixed(' '), '_'))
    for (m in names(cell_matrices)) {
      colnames(cell_matrices[[m]]) <- colnames(cell_matrices[[m]]) %>%
        stringr::str_remove(" \\(Opal .*\\)") %>%
        paste0("_", m)
    }
    cell_df <- do.call(cbind, c(list(cell_df), unname(cell_matrices)))

    # FIXME: Implement the addition of sample and tissue class variables more generally
    sample_df <- nc_read_data_frame(ds, 'sample')
    cell_df$t_number <- str_extract(fn, "T[0-9]+-[0-9]+(_[I0-9]+)?")
    cell_df$batch <- str_extract(fn, "batch[0-9]+")
    cell_df$panel <- sample_df$panel
    classifier_area <- nc_read_matrix(ds, 'area')
    stopifnot(nrow(classifier_area) == 1) # Only one sample
    cell_df$classifier_area <- classifier_area[1, as.character(cell_df$classifier_label)]

    # Cannot use on.exit, because it doesn't trigger on time with map_dfr
    if (is.character(fn)) {
      ncdf4::nc_close(ds)
    }

    cell_df
}

if (sys.nframe() == 0) {
  # Read all netCDF files on the command line and concatenate them into one big data frame
  # saved in an RDS file.
  parse_args <- function(args) {
    args <- commandArgs(T)
    stopifnot(length(args) >= 3)
    list(
      in_fns = args[1:(length(args) - 1)],
      out_fn = args[length(args)]
    )
  }

  args <- parse_args(commandArgs(T))

  objects <- rlang::set_names(args$in_fns) %>%
    purrr::map_dfr(~ read_cells(., intensity=T)) %>%
    dplyr::relocate(t_number, panel, .before = everything())

  saveRDS(objects, args$out_fn)
}
