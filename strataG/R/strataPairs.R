#' @name strataPairs
#' @title Strata Pairs
#' @description Make a matrix of all pairs of strata.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom utils combn
#' 
.strataPairs <- function(g) {
  if(nlevels(strata(g)) < 2) return(NULL)
  strata.vec <- sort(unique(as.character(strata(g))))
  strata.pairs <- t(combn(strata.vec, 2))
  colnames(strata.pairs) <- c("strata.1", "strata.2")
  as.data.frame(strata.pairs, stringsAsFactors = FALSE)
}