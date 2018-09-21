#' @title Number of Individuals Genotyped
#' @description Return the number of individuals genotyped for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results grouped by strata?
#'
#' @return vector of number of alleles per locus.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(msats.g)
#' 
#' numGenotyped(msats.g)
#'
#' @export
#' 
numGenotyped <- function(g, by.strata = FALSE) {
  df <- numMissing(g, by.strata = by.strata) 
  df <- if(by.strata) {
    df %>% 
      dplyr::left_join(
        g@data %>% 
          dplyr::group_by(stratum, locus) %>% 
          dplyr::summarize(num.ind = dplyr::n_distinct(id)),
        by = c("stratum", "locus")
      ) 
  } else {
    df %>% 
      dplyr::left_join(
        g@data %>% 
          dplyr::group_by(locus) %>% 
          dplyr::summarize(num.ind = dplyr::n_distinct(id)),
        by = c("locus")
      )
  }
  df %>% 
    dplyr::mutate(num.genotyped = num.ind - num.missing) %>% 
    dplyr::select(-num.missing, -num.ind) %>% 
    dplyr::ungroup() %>% 
    as.data.frame()
}