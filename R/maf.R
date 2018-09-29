#' @title Minimum Allele Frequencies
#' @description Calculate minimum allele frequencies for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results grouped by strata?
#'
#' @return A vector or matrix of minimum allele frequencies at each locus.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @seealso \link{alleleFreqs}
#'
#' @examples
#' data(msats.g)
#' 
#' maf(msats.g)
#' 
#' @export
#' 
maf <- function(g, by.strata = FALSE) {
  .calcMAF <- function(x) {
    maf <- min(x, na.rm = TRUE)
    if(maf == 1) 0 else maf
  }
    
  result <- alleleFreqs(g, by.strata = by.strata, type = "prop") %>% 
    sapply(function(loc) {
      if(by.strata) apply(loc, 2, .calcMAF) else .calcMAF(loc)
    })
  if(is.matrix(result)) t(result) else result
}