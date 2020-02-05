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
  cols <- if(by.strata) c("locus", "stratum") else "locus"
  as.data.frame(
    g@data[, 
           list(num.alleles = dplyr::n_distinct(allele, na.rm = TRUE)), 
           by = cols
          ]
  )
}