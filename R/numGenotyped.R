#' @title Number of Individuals Genotyped
#' @description Return the number of individuals genotyped for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results grouped by strata?
#' @param prop logical determining whether to return proportion missing.
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
numGenotyped <- function(g, by.strata = FALSE, prop = FALSE) {
  cols <- if(by.strata) c("locus", "stratum") else "locus"
  as.data.frame(
    if(prop) {
      g@data[, 
             list(num.genotyped = mean(!is.na(allele)) / getPloidy(g)), 
             by = cols
             ]
    } else {
      g@data[, 
             list(num.genotyped = sum(!is.na(allele)) / getPloidy(g)), 
             by = cols
             ]
    }
  )
}