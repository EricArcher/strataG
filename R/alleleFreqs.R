#' @title Allele Frequencies
#' @description Calculate allele frequencies for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results grouped by strata?
#'
#' @return A list of allele frequencies for each locus. Each element is a
#'   matrix or array with frequencies by count (\code{freq}) and 
#'   proportion (\code{prop}) of each allele.
#'   
#' @note If \code{g} is a haploid object with sequences, make sure to run 
#'   \code{\link{labelHaplotypes}} if sequences aren't already grouped by 
#'   haplotype.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @seealso \link{alleleFreqFormat}
#'
#' @examples
#' data(msats.g)
#' 
#' f <- alleleFreqs(msats.g)
#' f$D11t # Frequencies and proportions for Locus D11t
#' 
#' f.pop <- alleleFreqs(msats.g, TRUE)
#' f.pop$EV94[, "freq", "Coastal"] # Frequencies for EV94 in the Coastal population
#' 
#' @export
#' 
alleleFreqs <- function(g, by.strata = FALSE) {
  af <- if(by.strata) {
    g@data %>% 
      dplyr::group_by(stratum, locus, allele) %>% 
      dplyr::summarize(freq = n()) %>% 
      dplyr::filter(!is.na(allele)) %>% 
      dplyr::group_by(stratum, locus) 
  } else {
    g@data %>% 
      dplyr::group_by(locus, allele) %>% 
      dplyr::summarize(freq = n()) %>% 
      dplyr::filter(!is.na(allele)) %>% 
      dplyr::group_by(locus) 
  }
  af %>% 
    dplyr::mutate(prop = freq / sum(freq)) %>% 
    dplyr::ungroup()
}
