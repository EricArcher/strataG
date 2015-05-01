#' @title Expected Number of Alleles
#' @description Calculate the expected number of alleles in a sample of a
#'   given size and value of theta.
#'
#' @param n size of sample.
#' @param theta value of theta in population (= 2 * \code{ploidy} * Ne * mu).
#' @param ploidy ploidy of locus.
#'
#' @return a two element vector with the expected number of alleles
#'   (\code{num.alleles}) and variance (\code{var.num.alleles}).
#'
#' @references Ewens, W. 1972. The sampling theory of selectively neutral
#'   alleles. Theoretical Population Biology 3:87-112. Eqns. 11 and 24.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#'
#' exptdNumAlleles(20, 1, 2)
#'
#' # double the samples
#' exptdNumAlleles(40, 1, 2)
#'
#' # for a haploid locus
#' exptdNumAlleles(40, 1, 1)
#'
#' # double theta
#' exptdNumAlleles(40, 2, 1)
#'
#' @export
#'
exptdNumAlleles <- function(n, theta, ploidy) {
  n <- trunc(n)
  if(n < 1) stop("'n' must be greater than 1")
  if(theta <= 0) stop("'theta' must be greater than 0")
  ploidy <- trunc(ploidy)
  if(ploidy < 1) stop("'ploidy must be greater than 1")

  denom <- theta + 1:(ploidy * n - 1)
  e.k <- 1 + sum(theta / denom)
  var.e.k <- e.k - 1 - sum(theta ^ 2 / denom ^ 2)
  c(num.alleles = e.k, var.num.alleles = var.e.k)
}
