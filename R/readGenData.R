#' @title Read Genetic Data
#' @description A wrapper for \code{\link[data.table]{fread}} that sets 
#'   common values for missing data and removes blank lines.
#'
#' @param file filename of .csv file.
#' @param na.strings see \code{\link[data.table]{fread}}.
#' @param ... other arguments passed to \code{\link[data.table]{fread}}.
#'
#' @return a \code{data.frame}.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#' 
readGenData <- function(file, 
                        na.strings = c(NA, "NA", "", " ", "?", "."), ...) {
  df <- data.table::fread(
    file = file, 
    header = TRUE,
    na.strings = na.strings,
    colClasses = "character",
    ...
  ) 
  all.nas <- is.na(df)
  all.empty <- df == ""
  dplyr::filter(df, rowSums(all.nas | all.empty) != ncol(df))
}
