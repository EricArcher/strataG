#' @title Expected Heterozygosity
#' @description Calculate expected heterozygosity for diploid data.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{obsvdHet}}
#' 
#' @examples
#' data(dolph.msats)
#' data(dolph.strata)
#' msats.merge <- merge(dolph.strata[, c("ids", "fine")], dolph.msats, 
#'   all.y = TRUE)
#' msats.g <- df2gtypes(msats.merge, ploidy = 2)
#' exptdHet(msats.g)
#' 
#' @importFrom swfscMisc diversity
#' @export
#' 
exptdHet <- function(g) {
  apply(g@loci, 2, diversity)
}
