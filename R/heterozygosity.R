#' @title Observed and Expected Heterozygosity 
#' @description Calculate observed heterozygosity for diploid data.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results by strata?
#' @param type return \code{expected} or \code{observed} heterozygosity
#' 
#' @note For a measure of haplotypic diversity (haploid "heterozygosity"), 
#'   use \code{exptdHet}. If \code{g} is a haploid object with sequences, make sure to run 
#'   \code{\link{labelHaplotypes}} if sequences aren't already grouped by 
#'   haplotype.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' 
#' # Expected heterozygosity
#' exptdHet(msats.g)
#' 
#' # Observed heterozygosity
#' obsvdHet(msats.g)
#' 
heterozygosity <- function(g, by.strata = FALSE, type = c("expected", "observed")) {
  switch(
    match.arg(type),
    expected = .applyPerLocus(swfscMisc::diversity, g, by.strata = by.strata) %>% 
      dplyr::rename(exptd.het = value) %>% 
      as.data.frame(),
    observed = {
      is.het <- if(by.strata) {
        g@data %>% 
          dplyr::group_by(stratum, locus, id) %>% 
          dplyr::summarize(is.het = n_distinct(allele) > 1) %>% 
          dplyr::group_by(stratum, locus)
      } else {        
        g@data %>% 
          dplyr::group_by(locus, id) %>% 
          dplyr::summarize(is.het = n_distinct(allele) > 1) %>% 
          dplyr::group_by(locus)
      }
      is.het %>% 
        dplyr::summarize(obsvd.het = mean(is.het, na.rm = TRUE)) %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
    }
  )
}