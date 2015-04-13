#' @title Number Missing Data
#' @description Calculate the number of individuals with missing data by locus.
#'
#' @param g a \code{\link{gtypes}} object.
#' @param prop logical determining whether to return proportion missing.
#'
#' @return a vector of loci with number (or, if \code{prop = TRUE},
#'   the proportion) of individuals missing data for at least one allele.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' #msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' #msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5)
#' #num.missing(msats)
#' #num.missing(msats, prop = TRUE)
#'
#' @export

numMissing <- function(g, prop = FALSE) {
  apply(g@loci, 2, function(locus) {
    count <- sum(is.na(locus)) / g@ploidy
    if(prop) count <- count / nInd(g)
  })
}