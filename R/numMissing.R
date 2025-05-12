#' @title Number Missing Data
#' @description Calculate the number of individuals with missing data by locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results grouped by strata?
#' @param prop logical determining whether to return proportion missing.
#'
#' @return a vector of loci with number (or, if \code{prop = TRUE},
#'   the proportion) of individuals missing data for at least one allele.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(msats.g)
#' 
#' numMissing(msats.g)
#' numMissing(msats.g, prop = TRUE)
#'
#' @export
#' 
numMissing <- function(g, by.strata = FALSE, prop = FALSE) {
  cols <- if(by.strata) c("locus", "stratum") else "locus"
  as.data.frame(
    if(prop) {
      g@data[, 
             list(num.missing = mean(is.na(allele)) / getPloidy(g)), 
             by = cols
             ]
    } else {
      g@data[, 
             list(num.missing = sum(is.na(allele)) / getPloidy(g)), 
             by = cols
             ]
    }
  )
}