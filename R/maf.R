#' @title Minor Allele Frequencies
#' @description Calculate minor allele frequencies for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results grouped by strata?
#' @param maf.within if \code{by.strata = TRUE}, identify minor allele 
#'   within each strata independently? If \code{FALSE} minor allele is 
#'   identified from all individuals.
#'
#' @return A vector or matrix of minor allele frequencies at each locus.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @seealso \link{alleleFreqs}
#'
#' @examples
#' data(msats.g)
#' 
#' maf(msats.g)
#' 
#' # minor allele identified from all indivudals
#' maf(msats.g, by.strata = TRUE)
#' 
#' # minor allele identified within each strata
#' maf(msats.g, by.strata = TRUE, maf.within = TRUE)
#' 
#' @export
#' 
maf <- function(g, by.strata = FALSE, maf.within = FALSE) {
  .calcMAF <- function(x) {
    maf <- min(x, na.rm = TRUE)
    if(maf == 1) 0 else maf
  }
  
  if(by.strata) {
    result <- if(maf.within) {
      g %>% 
        alleleFreqs(by.strata = TRUE, type = "prop") %>% 
        purrr::map(function(x) apply(x, 2, .calcMAF))
    } else {
      st.af <- alleleFreqs(g, by.strata = TRUE, type = "prop")
      alleleFreqs(g, by.strata = FALSE, type = "prop") %>% 
        purrr::imap(function(freqs, locus) {
          allele <- names(freqs)[which.min(freqs)]
          st.af[[locus]][allele, ]
        })
    }
    do.call(rbind, result) %>% 
      as.data.frame %>% 
      tibble::rownames_to_column("locus") %>% 
      tidyr::gather("stratum", "maf", -.data$locus) %>% 
      dplyr::select(.data$stratum, .data$locus, .data$maf)
  } else {
    alleleFreqs(g, by.strata = FALSE, type = "prop") %>% 
      purrr::map_dbl(.calcMAF) %>%  
      utils::stack() %>% 
      stats::setNames(c("maf", "locus")) %>% 
      dplyr::select(.data$locus, .data$maf)
  }
}