#' @title Read Genetic Data
#' @description A wrapper for \code{\link{read.csv}} that sets common values
#'   for missing data and removes blank lines.
#'
#' @param file filename of .csv file.
#' @param na.strings see \code{\link{read.table}}.
#' @param ... other arguments passed to \code{\link{read.table}}.
#'
#' @return a \code{data.frame}.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#' 
readGenData <- function(file, na.strings = c(NA, "NA", "", " ", "?", "."), ...) {
  utils::read.csv(
    file = file, 
    na.strings = na.strings,
    colClasses = "character", 
    stringsAsFactors = FALSE, 
    ...
  ) %>% 
    dplyr::filter(rowSums(is.na(.data) | .data == "") != ncol(.data))
}
