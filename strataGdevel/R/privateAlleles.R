#' @title Private Alleles
#' @description The number of private alleles in each strata and locus.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return matrix with the number of private alleles in each strata at each locus.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.msats)
#' data(dolph.strata)
#' msats.merge <- merge(dolph.strata[, c("ids", "fine")], dolph.msats, all.y = TRUE)
#' msats <- df2gtypes(msats.merge, ploidy = 2)
#' 
#' privateAlleles(msats)
#' 
#' @export
#' 
privateAlleles <- function(g) {
  freqs <- alleleFreqs(g, T)
  
  do.call(rbind, lapply(freqs, function(f) {
    f <- f[, 1, ]
    f[f > 0] <- 1
    pa <- rowSums(apply(f, 1, function(x) {
      if(sum(x > 0) == 1) x else rep(0, length(x))
    }))
    names(pa) <- colnames(f)
    pa
  }))
}