#' @title Nm from Fst
#' @description Calculate Nm (number of migrants per generation) for a 
#'   given value of Fst.
#'
#' @param fst estimate of Fst between populations.
#' @param ploidy ploidy of locus Fst is from.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @seealso \link{wrightFst}, \link{numGensEq}, \link{expectedNumAlleles}
#'
#' @examples 
#' fst <- seq(0.001, 0.2, length.out = 100)
#' Nm <- fstToNm(fst, 2)
#' plot(fst, Nm, type = "l")
#'
#' @export
#' 
fstToNm <- function(fst, ploidy) {
  ((1 / fst) - 1) / (ploidy * 2)
}