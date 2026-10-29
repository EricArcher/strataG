#' @title Read Genetic Data
#' @description A wrapper for \code{\link[data.table]{fread}} that sets 
#'   common values for missing data and removes blank lines.
#'
#' @param file filename of .csv file.
#' @param na.strings a vector of strings that are interpreted as missing values.
#'   See \code{\link[data.table]{fread}} for more information.
#' @param ... other arguments passed to \code{\link[data.table]{fread}}.
#'
#' @return a \code{data.frame}.
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#'
#' @export
#' 
readGenData <- function(
 file, na.strings = c("NA", "", "?", "."), ...
) {
  df <- data.table::fread(
    file = file, 
    header = TRUE,
    strip.white = TRUE,
    colClasses = "character",
    na.strings = na.strings,
    ...
  )
  
  all.nas <- is.na(df)
  all.empty <- df == ""
  
  df |> 
    dplyr::filter(rowSums(all.nas | all.empty) != ncol(df)) |> 
    as.data.frame()
    
}
