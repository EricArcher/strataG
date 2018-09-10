#' @title Proportion Unique Alleles
#' @description Calculate the proportion of alleles that are unique.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results by strata?
#' 
#' @return a vector of the proportion of unique (occuring only in one individual) 
#'   alleles for each locus.
#' 
#' @note If \code{g} is a haploid object with sequences, make sure to run 
#'   \code{\link{labelHaplotypes}} if sequences aren't already grouped by 
#'   haplotype.
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
propUniqueAlleles <- function(g, by.strata = TRUE) { 
  df <- if(by.strata) {
    g@data %>% 
      group_by(stratum, locus, allele) %>% 
      summarize(n = n_distinct(id)) %>% 
      group_by(stratum, locus) %>% 
      summarize(num.unique = sum(n == 1)) %>% 
      left_join(numGenotyped(g, by.strata), by = c("stratum", "locus")) 
  } else {
    g@data %>% 
      group_by(locus, allele) %>% 
      summarize(n = n_distinct(id)) %>% 
      group_by(locus) %>% 
      summarize(num.unique = sum(n == 1)) %>% 
      left_join(numGenotyped(g, by.strata), by = c("locus"))
  } 
  df %>% 
    mutate(prop.unique.alleles = num.unique / num.genotyped) %>% 
    select(-num.unique, -num.genotyped) %>% 
    ungroup
}