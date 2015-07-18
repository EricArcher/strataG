#' @name heterozygosity
#' @title Observed and Expected Heterozygosity 
#' @description Calculate observed heterozygosity for diploid data.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.msats)
#' data(dolph.strata)
#' msats.merge <- merge(dolph.strata[, c("ids", "fine")], dolph.msats,  all.y = TRUE)
#' msats <- df2gtypes(msats.merge, ploidy = 2)
#' 
#' obsvdHet(msats)
#' exptdHet(msats)
#' 
#' @export
#' 
obsvdHet <- function(g) {
  apply(loci(g), 2, function(locus) {
    locus.mat <- na.omit(matrix(locus, ncol = ploidy(g)))
    is.homozgt <- apply(locus.mat, 1, function(x) length(unique(x)) == 1)
    1 - (sum(is.homozgt) / nrow(locus.mat))
  })
}


#' @rdname heterozygosity
#' @importFrom swfscMisc diversity
#' @export
#' 
exptdHet <- function(g) {
  apply(loci(g), 2, swfscMisc::diversity)
}