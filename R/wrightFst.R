#' @title Wright's Fst
#' @description Calcualte Wright's Fst from Ne, dispersal, and generation time.
#' 
#' @param Ne Effective population size.
#' @param dispersal migration rate in terms of probability of an individual 
#'   migrating in a generation.
#' @param gen.time number of generations since ancestral population.
#' @param ploidy ploidy of the locus
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{numGensEq}, \link{fstToNm}, \link{expectedNumAlleles}
#' 
#' @examples 
#' dispersal <- seq(0.05, 0.8, by = 0.05)
#' fst <- wrightFst(100, dispersal, 20, 2)
#' plot(dispersal, fst, type = "l")
#' 
#' @export
#' 
wrightFst <- function(Ne, dispersal, gen.time, ploidy) {
  1 / (2 * ploidy * Ne * dispersal * gen.time + 1)
}