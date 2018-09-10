#' @title Allelic Richness
#' @description Calculate allelic richness for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results grouped by strata?
#' 
#' @return the allelic richness of each locus calculated as the number of 
#'   alleles divided by the number of samples without missing data at 
#'   that locus.
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
  if(ploidy(g) == 1 & !is.null(sequences(g))) g <- labelHaplotypes(g)$gtypes
  join.by <- if(by.strata) c("stratum", "locus") else "locus"
  numAlleles(g, by.strata) %>% 
    dplyr::left_join(numGenotyped(g, by.strata), by = join.by) %>% 
    dplyr::mutate(allelic.richness = num.alleles / num.genotyped) %>% 
    dplyr::select(-num.alleles, -num.genotyped) %>% 
    as.data.frame()
}
