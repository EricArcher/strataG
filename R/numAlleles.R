#' @title Number of Alleles
#' @description Return the number of alleles for each locus.
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
#' numAlleles(msats.g)
#'
#' @export
#' 
numAlleles <- function(g, by.strata = FALSE) {
  g <- .checkHapsLabelled(g)
  cols <- if(by.strata) c("stratum", "locus") else "locus"
  result <- g@data[, 
                   list(num.alleles = dplyr::n_distinct(allele, na.rm = TRUE)), 
                   by = cols]
  result <- if(by.strata) {
    dplyr::arrange(result, .data$stratum, .data$locus)
  } else {
    dplyr::arrange(result, .data$locus)
  }
  as.data.frame(result)
}