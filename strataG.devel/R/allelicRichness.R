#' @title Allelic Richness
#' @description Calculate allelic richness for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export
#' 
allelicRichness <- function(g) {
  apply(g@loci, 2, function(locus) {
    locus <- na.omit(locus)
    length(unique(locus)) / (length(locus) / g@ploidy)
  })
}
