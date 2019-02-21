#' @title Minimum Allele Frequencies
#' @description Calculate minimum allele frequencies for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results grouped by strata?
#'
#' @return A vector or matrix of minimum allele frequencies at each locus.
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
#' @export
#' 
maf <- function(g, by.strata = FALSE) {
  .calcMAF <- function(x) {
    maf <- min(x, na.rm = TRUE)
    if(maf == 1) 0 else maf
  }
  
  if(by.strata) {
    result <- g %>% 
      alleleFreqs(by.strata = TRUE, type = "prop") %>% 
      purrr::map(function(x) apply(x, 2, .calcMAF))
    result <- do.call(rbind, result) %>% 
      as.data.frame %>% 
      tibble::rownames_to_column("locus") %>% 
      tidyr::gather("stratum", "maf", -.data$locus) %>% 
      dplyr::select(.data$stratum, .data$locus, .data$maf)
  } else {
    g %>% 
      alleleFreqs(by.strata = FALSE, type = "prop") %>% 
      purrr::map_dbl(.calcMAF) %>%  
      utils::stack() %>% 
      stats::setNames(c("maf", "locus")) %>% 
      dplyr::select(.data$locus, .data$maf)
  }
}