#' @title Allele Frequencies
#' @description Calculate allele frequencies for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical. If \code{TRUE} every element in the return list is 
#'   a three dimensional array where the third dimension contains frequencies 
#'   and proportions for each stratum.
#'
#' @return A list of allele frequencies for each locus. Each element is a
#'   matrix or array with frequencies by count (\code{freq}) and 
#'   proportion (\code{prop}) of each allele.
#'   
#' @note If \code{g} is a haploid object with sequences, make sure to run 
#'   \code{\link{labelHaplotypes}} if sequences aren't already grouped by 
#'   haplotype.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @seealso \link{alleleFreqFormat}
#'
#' @examples
#' data(msats.g)
#' 
#' f <- alleleFreqs(msats.g)
#' f$D11t # Frequencies and proportions for Locus D11t
#' 
#' f.pop <- alleleFreqs(msats.g, TRUE)
#' f.pop$EV94[, "freq", "Coastal"] # Frequencies for EV94 in the Coastal population
#' 
#' @export

alleleFreqs <- function(g, by.strata = FALSE) {
  freqs <- vector("list", length = nLoc(g))
  names(freqs) <- locNames(g)
  if(by.strata & nStrata(g) > 1) {
    for(i in locNames(g)) {
      f <- table(g@data[[i]], g@data$strata)
      p <- prop.table(f, 2)
      freqs[[i]] <- array(
        dim = list(nrow(f), 2, ncol(f)),
        dimnames = list(rownames(f), c("freq", "prop"), colnames(f))
      )
      for(j in 1:ncol(f)) freqs[[i]][, , j] <- cbind(f[, j], p[, j])
    } 
  } else {
    for(i in locNames(g)) {
      f <- table(g@data[[i]])
      freqs[[i]] <- cbind(freq = f, prop = f / sum(f))
    }
  }
  freqs
}
