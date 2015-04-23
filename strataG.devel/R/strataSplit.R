#' @title Split Strata
#' @description Return a list of gtypes for each strata
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param strata a character vector giving a subset of strata to select. 
#'   If \code{NULL} then a list with all strata is created.
#' 
#' @return A named list where each element is a \code{gtypes} object 
#'   for a single stratum in \code{g}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
strataSplit <- function(g, strata = NULL) {
  if(is.null(strata)) strata <- strataNames(g)
  g <- subset(g, strata)
  sapply(strata, function(st) subset(g, st), simplify = FALSE)
}