#' @title Number of Alleles
#' @description Return the number of alleles for each locus.
#'
#' @param g a \code{\link{gtypes}} object.
#'
#' @return vector of number of alleles per locus.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(dolph.msats)
#' data(dolph.strata)
#' msats.merge <- merge(dolph.strata[, c("ids", "fine")], dolph.msats, all.y = TRUE)
#' msats <- df2gtypes(msats.merge, ploidy = 2)
#' 
#' numAlleles(msats)
#'
#' @export
#' 
numAlleles <- function(g) {
  apply(g@loci, 2, function(locus) length(unique(na.omit(locus))))
}