#' @title Allele Frequencies
#' @description Calculate allele frequencies or proportions for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results by strata?
#' @param type return counts (\code{"freq"}) or proportions (\code{"prop"})
#'
#' @return A list of allele frequencies for each locus. Each element is a
#'   vector (\code{by.strata = FALSE}) or matrix (\code{by.strata = TRUE}) 
#'   with the frequency or proportion of each allele. 
#'   
#' @note If \code{g} is a haploid object with sequences, make sure to run 
#'   \code{\link{labelHaplotypes}} if sequences aren't already grouped by 
#'   haplotype.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(msats.g)
#' 
#' f <- alleleFreqs(msats.g)
#' f$D11t # Frequencies for Locus D11t
#' 
#' f.pop <- alleleFreqs(msats.g, TRUE, "prop")
#' f.pop$EV94[, "Coastal"] # Proportions of EV94 alleles in the Coastal population
#' 
#' @export
#' 
alleleFreqs <- function(g, by.strata = FALSE, type = c("freq", "prop")) {
  af <- if(by.strata) {
    split(g@data, g@data$locus) %>% 
      purrr::map(function(x) table(x$allele, x$stratum))
  } else {
    split(g@data, g@data$locus) %>% 
      purrr::map(function(x) table(x$allele))
  }
  
  if(match.arg(type) == "prop") {
    af %>% 
      purrr::map(function(x) {
        if(length(dim(x)) == 1) {
          prop.table(x)
        } else {
          prop.table(x, length(dim(x)))
        }
      })
  } else {
    af
  }
}
