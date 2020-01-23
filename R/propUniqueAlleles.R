#' @title Proportion Unique Alleles
#' @description Calculate the proportion of alleles that are unique.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results grouped by strata?
#' 
#' @return a vector of the proportion of unique (occuring only in one individual) 
#'   alleles for each locus.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{privateAlleles}
#' 
#' @examples
#' data(msats.g)
#' 
#' propUniqueAlleles(msats.g)
#' 
#' @export
#' 
propUniqueAlleles <- function(g, by.strata = FALSE) { 
  g <- .checkHapsLabelled(g)
  
  if(by.strata) {
    g@data %>% 
      dplyr::group_by(.data$stratum, .data$locus, .data$allele) %>% 
      dplyr::summarize(n = dplyr::n_distinct(.data$id)) %>% 
      dplyr::ungroup() %>% 
      dplyr::group_by(.data$stratum, .data$locus) %>% 
      dplyr::summarize(num.unique = sum(.data$n == 1)) %>% 
      dplyr::ungroup() %>% 
      dplyr::left_join(numGenotyped(g, by.strata), by = c("stratum", "locus")) 
  } else {
    g@data %>% 
      dplyr::group_by(.data$locus, .data$allele) %>% 
      dplyr::summarize(n = dplyr::n_distinct(.data$id)) %>% 
      dplyr::ungroup() %>% 
      dplyr::group_by(.data$locus) %>% 
      dplyr::summarize(num.unique = sum(.data$n == 1)) %>% 
      dplyr::ungroup() %>% 
      dplyr::left_join(numGenotyped(g, by.strata), by = c("locus"))
  } %>% 
    dplyr::mutate(prop.unique.alleles = .data$num.unique / .data$num.genotyped) %>% 
    dplyr::select(-.data$num.unique, -.data$num.genotyped) %>% 
    as.data.frame()
}