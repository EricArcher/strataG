#' @title CHI-squared Analysis of Population Structure
#'
#' @param g a \code{gtypes} object.
#'
#' @return a \code{\link{gtype.struct.stat}} list with the CHI-squared estimate.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(dolph.strata)
#' data(dolph.msats)
#' #msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE, sort = FALSE)
#' #msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "msats")
#'
#' #stat.chi2(msats)
#'
#' @export

statChi2 <- function(g) {
  chi2 <- vector("numeric", ncol(g@locus.data))
  strata <- rep(g@strata, g@ploidy)
  for(i in 1:length(chi2)) {
    obs.freq <- table(g@locus.data[, i], strata, useNA = "no")
    if(nrow(obs.freq) > 0 & ncol(obs.freq) > 1) {
      exp.freq <- outer(rowSums(obs.freq), colSums(obs.freq)) / sum(obs.freq)
      chi2[i] <- sum((obs.freq - exp.freq) ^ 2 / exp.freq)
    }
  }
  est <- sum(chi2, na.rm = TRUE)

  list(stat.name = "Chi2", estimate = est, 
       strata.freq = table(strata, useNA = "no"))
}