#' @title Locus Summaries
#' @description Compile standard by-locus summaries.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical. If \code{TRUE}, return a list of summary matrices 
#'   for each stratum.
#' 
#' @return A matrix with rows for each locus and columns containing summaries of:
#' \describe{
#'   \item{\code{num.genotyped}}{The number of samples genotyped}
#'   \item{\code{prop.genotyped}}{The proportion of samples genotyped}
#'   \item{\code{num.alleles}}{The number of alleles in the locus}
#'   \item{\code{allelic.richness}}{The allelic richness of the locus}
#'   \item{\code{prop.unique.alleles}}{Proportion of alleles found in a single sample}
#'   \item{\code{expt.heterozygosity}}{Expected heterozygosity}
#'   \item{\code{obsvd.heterozygosity}}{Observed heterozygosity}
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#' 
#' summarizeLoci(msats.g)
#' 
#' @export
#' 
summarizeLoci <- function(g, by.strata = FALSE) {
  by.cols <- if(by.strata) c("stratum", "locus") else "locus"
  smry <- numGenotyped(g, by.strata) %>% 
    dplyr::left_join(numMissing(g, by.strata), by = by.cols) %>% 
    dplyr::mutate(prop.genotyped = .data$num.genotyped / (.data$num.genotyped + .data$num.missing)) %>% 
    dplyr::left_join(numAlleles(g, by.strata), by = by.cols) %>% 
    dplyr::left_join(allelicRichness(g, by.strata), by = by.cols) %>% 
    dplyr::left_join(propUniqueAlleles(g, by.strata), by = by.cols) %>% 
    dplyr::left_join(heterozygosity(g, by.strata, "expected"), by = by.cols) %>% 
    dplyr::left_join(heterozygosity(g, by.strata, "observed"), by = by.cols) 

  if(getPloidy(g) == 1) {
    smry <- smry %>% 
      dplyr::rename(
        num.haplotypes = .data$num.alleles,
        prop.unique.haplotypes = .data$prop.unique.alleles,
        haplotypic.diversity = .data$exptd.het.x
      ) %>% 
      dplyr::select(-.data$exptd.het.y)
  }
  
  smry
}