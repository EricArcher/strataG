#' @title Allelic Richness
#' @description Calculate allelic richness for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results grouped by strata?
#' 
#' @return a data frame with the allelic richness of each locus calculated
#'   as the number of alleles divided by the number of samples without
#'   missing data at that locus.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples 
#' data(msats.g)
#' allelicRichness(msats.g)
#'
#' @export
#' 
allelicRichness <- function(g, by.strata = FALSE) {
  g <- .checkHapsLabelled(g)
  join.by <- if(by.strata) c("stratum", "locus") else "locus"
  # average number of alleles per locus genotyped
  numAlleles(g, by.strata) %>% 
    dplyr::left_join(numGenotyped(g, by.strata), by = join.by) %>% 
    dplyr::mutate(allelic.richness = .data$num.alleles / .data$num.genotyped) %>% 
    dplyr::select(-.data$num.alleles, -.data$num.genotyped) %>% 
    as.data.frame()
}