#' Dichotomize column around the median (greater than or equal to)
#'
#' @param df A data frame
#' @param from.col Character, the column name containing a numeric variable
#' @param to.col Character, the name of the column containing the median category
#' @param verbose Whether to print the column median
#'
#' @return A data frame with a new column containing a factor with levels "high"
#' (greater than or equal to the median) or "low" (less than the median)
#' @export
#'
#' @examples dichotomize_median(patientdata, from.col = "TI_CD38_total", to.col = "TI_CD38"
dichotomize_median <- function(df, from.col, to.col, verbose = T){
  
  med <- median(df[[from.col]], na.rm=T)
  
  if(verbose){print(paste("Median of", from.col, "is", med))}
  
  df[[to.col]] <- factor(ifelse(df[[from.col]] >= med, "high", "low"),
                         levels = c("low", "high"))
  df %>%
    relocate(all_of(to.col), .after = all_of(from.col))
}

#' Create quantile categories in a new column
#'
#' @param df A data frame
#' @param from.col Character, the column name containing a numeric variable
#' @param to.col Character, the name of the column containing the quantile category
#' @param n The number of quantiles. Must be 2, 3 or 4.
#'
#' @return A data frame with a new column containing a factor based on the new quantiles.
#' @export
#'
#' @examples 
create_ntiles <- function(df, from.col, to.col, n = 3){
  
  stopifnot(n %in% c(2,3,4))
  
  labels = switch(as.character(n),
                  "2" = c("low", "high"),
                  "3" = c("low", "medium", "high"),
                  "4" = c("Q1", "Q2", "Q3", "Q4")
  )
  
  df[[to.col]] <- factor(dplyr::ntile(df[[from.col]], n), labels = labels)
  
  df %>%
    relocate(all_of(to.col), .after = all_of(from.col))
}