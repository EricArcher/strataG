#' @title Allelic Richness
#' @description Calculate allelic richness for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return the allelic richness of each locus calculated as the number of 
#'   alleles divided by the number of samples without missing data at 
#'   that locus.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples 
#' data(dolph.msats)
#' data(dolph.strata)
#' msats.merge <- merge(dolph.strata[, c("ids", "fine")], dolph.msats, all.y = TRUE)
#' msats <- df2gtypes(msats.merge, ploidy = 2)
#' 
#' allelicRichness(msats)
#'
#' @export
#' 
allelicRichness <- function(g) {
  apply(g@loci, 2, function(locus) {
    locus <- na.omit(locus)
    length(unique(locus)) / (length(locus) / ploidy(g))
  })
}
