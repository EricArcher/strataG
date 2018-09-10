#' @title Number of Alleles
#' @description Return the number of alleles for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results by strata?
#'
#' @return vector of number of alleles per locus.
#'
#' @note If \code{g} is a haploid object with sequences, make sure to run 
#'   \code{\link{labelHaplotypes}} if sequences aren't already grouped by 
#'   haplotype.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(msats.g)
#' 
#' numAlleles(msats.g)
#'
#' @importFrom stats na.omit
#' @export
#' 
numAlleles <- function(g, by.strata = TRUE) {
  .countAlleles <- function(x) n_distinct(x, na.rm = TRUE)
  .applyPerLocus(.countAlleles, g, by.strata = by.strata) %>%
    rename(num.alleles = value)
}