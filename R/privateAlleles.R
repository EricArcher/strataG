#' @title Private Alleles
#' @description The number of private alleles in each strata and locus.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return matrix with the number of private alleles in each strata at each 
#'   locus. This is the number of alleles only present in one stratum.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{propUniqueAlleles}
#' 
#' @examples
#' data(msats.g)
#' 
#' privateAlleles(msats.g)
#' 
#' @export
#' 
privateAlleles <- function(g) {
  do.call(
    rbind,
    sapply(
      alleleFreqs(.checkHapsLabelled(g), by.strata = TRUE), 
      function(f) {
        f[f > 0] <- 1
        apply(f, 1, function(x) {
          if(sum(x > 0) == 1) x else rep(0, length(x))
        }) %>% 
          rbind() %>% 
          rowSums() %>% 
          stats::setNames(dimnames(f)[[2]])
      },
      simplify = FALSE
    )
  )
}