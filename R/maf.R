#' @title Minimum Allele Frequencies
#' @description Calculate minimum allele frequencies for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical. If \code{TRUE} every element in the return list is 
#'   a three dimensional array where the third dimension contains frequencies 
#'   and proportions for each stratum.
#'
#' @return A vector or array of minimum allele frequencies at each locus.
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

maf <- function(g, by.strata = FALSE) {
  freqs <- alleleFreqs(g, by.strata)
  sapply(freqs, function(loc) {
    if(by.strata) {
      apply(loc, 3, function(x) {
        maf <- min(x[, "prop"])
        if(maf == 1) 0 else maf
      })
    } else {
      maf <- min(loc[, "prop"])
      if(maf == 1) 0 else maf
    }
  })
}