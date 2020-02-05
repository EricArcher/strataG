#' @title Theta
#' @description Calculate theta from heterozygosity of each locus.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results grouped by strata?
#' 
#' @return vector of theta values for each locus.
#' 
#' @details Calculates theta for each locus using the
#'   \code{\link[pegas]{theta.h}} function.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' 
#' theta(msats.g)
#' 
#' @export
#' 
theta <- function(g, by.strata = FALSE) {
  .thetaFunc <- function(x) pegas::theta.h(stats::na.omit(x))
  .applyPerLocus(.thetaFunc, g, by.strata = by.strata) %>%
    dplyr::rename(theta = .data$value) %>% 
    as.data.frame()
}