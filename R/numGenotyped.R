#' @title Number of Individuals Genotyped
#' @description Return the number of individuals genotyped for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#'
#' @return vector of number of alleles per locus.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(msats.g)
#' 
#' numGenotyped(msats.g)
#'
#' @export
#' 
numGenotyped <- function(g) {
  nInd(g) - numMissing(g)
}