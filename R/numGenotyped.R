#' @title Number of Individuals Genotyped
#' @description Return the number of individuals genotyped for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results by strata?
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
numGenotyped <- function(g, by.strata = TRUE) {
  df <- numMissing(g, by.strata = by.strata) 
  df <- if(by.strata) {
    df %>% 
      left_join(
        g@data %>% 
          group_by(stratum, locus) %>% 
          summarize(num.ind = n_distinct(id)),
        by = c("stratum", "locus")
      ) 
  } else {
    df %>% 
      left_join(
        g@data %>% 
          group_by(locus) %>% 
          summarize(num.ind = n_distinct(id)),
        by = c("locus")
      )
  }
  df %>% 
    mutate(num.genotyped = num.ind - num.missing) %>% 
    select(-num.missing, -num.ind) %>% 
    ungroup
}