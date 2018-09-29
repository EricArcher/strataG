#' @title Population Genetics Equations
#' @description Collection of classical population genetics equations.
#' 
#' @param Ne Effective population size.
#' @param dispersal Migration rate in terms of probability of an individual 
#'   migrating in a generation.
#' @param fst value of Fst at equilibrium.
#' @param gen.time Number of generations since ancestral population.
#' @param ploidy Ploidy of the locus.
#' @param n Sample size.
#' @param theta Product of effective population size (Ne) and mutation rate (mu).
#' 
#' @details \describe{
#'   \item{wrightFst}{Calculate Wright's Fst from Ne, dispersal, and generation time.}
#'   \item{numGensEq}{Calculate the number of generations to equilibrium based on a an
#'   ideal Wright model.}
#'   \item{fstToNm}{Calculate Nm (number of migrants per generation) for a 
#'   given value of Fst.}
#'   \item{expectedNumAlleles}{Calculate the expected number of alleles in a sample of a
#'   given size and value of theta.}
#' }
#' 
#' @return \describe{
#'   \item{wrightFst}{}
#'   \item{numGensEq}{}
#'   \item{fstToNm}{}
#'   \item{expectedNumAlleles}{a two element vector with the expected number of alleles
#'   (\code{num.alleles}) and variance (\code{var.num.alleles}).}
#' }
#'
#' @references Ewens, W. 1972. The sampling theory of selectively neutral
#'   alleles. Theoretical Population Biology 3:87-112. Eqns. 11 and 24.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples 
#' dispersal <- seq(0.05, 0.8, by = 0.05)
#' fst <- wrightFst(100, dispersal, 20, 2)
#' plot(dispersal, fst, type = "l")
#' 
#' numGensEq(0.15, 100, 20)
#' numGensEq(0.3, 100, 20)
#' numGensEq(0.15, 50, 20)
#' 
#' fst <- seq(0.001, 0.2, length.out = 100)
#' Nm <- fstToNm(fst, 2)
#' plot(fst, Nm, type = "l")
#' 
#' expectedNumAlleles(20, 1, 2)
#' # double the samples
#' expectedNumAlleles(40, 1, 2)
#' # for a haploid locus
#' expectedNumAlleles(40, 1, 1)
#' # double theta
#' expectedNumAlleles(40, 2, 1)
#'
#' @name popGenEqns
NULL


#' @rdname popGenEqns
#' @export
#' 
wrightFst <- function(Ne, dispersal, gen.time, ploidy) {
  1 / (2 * ploidy * Ne * dispersal * gen.time + 1)
}


#' @rdname popGenEqns
#' @export
#' 
numGensEq <- function(fst, Ne, gen.time) {
  term1 <- log(1 - fst)
  term2 <- 1 - (1 / (2 * Ne))
  n.gens <- term1 / log(term2)
  eq.fst <- 1 - (term2 ^ gen.time)
  cbind(n.gens = n.gens, eq.fst = eq.fst)
}


#' @rdname popGenEqns
#' @export
#' 
fstToNm <- function(fst, ploidy) {
  ((1 / fst) - 1) / (ploidy * 2)
}


#' @rdname popGenEqns
#' @export
#' 
expectedNumAlleles <- function(n, theta, ploidy) {
  n <- trunc(n)
  ploidy <- trunc(ploidy)
  result <- c(num.alleles = NA, var.num.alleles = NA)
  
  if(n < 1) {
    warning("'n' must be 1 or greater. NA returned")
    return(result)
  }
  if(theta <= 0) {
    warning("'theta' must be greater than 0. NA returned")
    return(result)
  }
  if(ploidy < 1) {
    warning("'ploidy' must be 1 or greater. NA returned")
    return(result)
  }
  
  denom <- theta + 1:(ploidy * n - 1)
  result["num.alleles"] <- 1 + sum(theta / denom)
  var.term <- 1 - sum(theta ^ 2 / denom ^ 2)
  result["var.num.alleles"] <- result["num.alleles"] - var.term
  result
}