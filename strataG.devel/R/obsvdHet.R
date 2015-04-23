#' @title Observed Heterozygosity 
#' @description Calculate observed heterozygosity for diploid data.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{exptdHet}}
#' 
#' @examples
#' data(dolph.msats)
#' data(dolph.strata)
#' msats.merge <- merge(dolph.strata[, c("ids", "fine")], dolph.msats, 
#'   all.y = TRUE)
#' msats.g <- df2gtypes(msats.merge, ploidy = 2)
#' obsvdHet(msats.g)
#' 
#' @export
#' 
obsvdHet <- function(g) {
  apply(g@loci, 2, function(locus) {
    locus.mat <- na.omit(matrix(locus, ncol = g@ploidy))
    is.homozgt <- apply(locus.mat, 1, function(x) length(unique(x)) == 1)
    1 - (sum(is.homozgt) / nrow(locus.mat))
  })
}