library(ncdf4)
library(tibble)

pixel_size <- 0.0004286127  # mm^2  # FIXME: Check with Iris if this is corrrect
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

read_cells <- function(fn) {
    if (is.character(fn)) {
      ds <- ncdf4::nc_open(fn)
      on.exit(ncdf4::nc_close(ds), add=T)
    } else {
      ds <- fn
    }
    df <- nc_read_data_frame(ds, "cell")
    positive <- nc_read_matrix(ds, "positive_classification") > 0

    positive <- as_tibble(t(positive), rownames = "cell") %>%
        dplyr::rename_at(2:ncol(.), ~ paste0(stringr::str_extract(., "^[A-Za-z0-9]+"), "_positive")) %>%
        dplyr::mutate(cell = as.integer(cell))
    df <- df %>%
        dplyr::mutate(
            x = pixel_size * (xmin + xmax) / 2,
            y = pixel_size * (ymin + ymax) / 2) %>%
        dplyr::distinct()
    df <- dplyr::left_join(df, positive, by = "cell")

    df
}
