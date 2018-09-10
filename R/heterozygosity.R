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
heterozygosity <- function(g, by.strata = TRUE, type = c("expected", "observed")) {
  type = match.arg(type)
  het <- switch(
    type,
    expected = .applyPerLocus(swfscMisc::diversity, g, by.strata = by.strata) %>% 
      rename(exptd.het = value),
    observed = {
      is.het <- if(by.strata) {
        g@data %>% 
          group_by(stratum, locus, id) %>% 
          summarize(is.het = n_distinct(allele) > 1) %>% 
          group_by(stratum, locus)
      } else {        
        g@data %>% 
          group_by(locus, id) %>% 
          summarize(is.het = n_distinct(allele) > 1) %>% 
          group_by(locus)
      }
      is.het %>% 
        summarize(obsvd.het = mean(is.het, na.rm = TRUE)) %>% 
        ungroup
    }
  )
}