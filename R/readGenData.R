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
                        na.strings = c("NA", "", "?", "."), ...) {
  .replaceNA <- function(x, na.strings) {
    ifelse(x %in% na.strings, NA, x)
  }
  
  df <- data.table::fread(
    file = file, 
    header = TRUE,
    strip.white = TRUE,
    colClasses = "character",
    ...
  ) %>% 
    as.data.frame() %>% 
    dplyr::mutate(dplyr::across(
      .cols = dplyr::everything(), 
      .fns = .replaceNA, 
      na.strings = na.strings
    ))
  all.nas <- is.na(df)
  all.empty <- df == ""
  dplyr::filter(df, rowSums(all.nas | all.empty) != ncol(df))
}