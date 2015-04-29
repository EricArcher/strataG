#' @title Strata Pairs
#' @description Make a matrix of all pairs of strata.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
strataPairs <- function(g) {
  st <- levels(strata(g))
  if(length(st) < 2) stop("'g' has fewer than 2 strata")
  mat <- t(combn(strata, 2))
  colnames(mat) <- c("strata.1", "strata.2")
  mat
}