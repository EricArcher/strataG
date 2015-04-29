#' @rdname popStructStat
#' @export
#' 
statChi2 <- function(g, strata = NULL, ...) {
  strata <- if(is.null(strata)) {
    rep(g@strata, g@ploidy)
  } else {
    rep(rep(strata, length.out = nInd(g)), g@ploidy)
  }
  
  chi2 <- vector("numeric", ncol(g@loci))
  for(i in 1:length(chi2)) {
    obs.freq <- table(g@loci[, i], strata, useNA = "no")
    if(nrow(obs.freq) > 0 & ncol(obs.freq) > 1) {
      exp.freq <- outer(rowSums(obs.freq), colSums(obs.freq)) / sum(obs.freq)
      chi2[i] <- sum((obs.freq - exp.freq) ^ 2 / exp.freq, na.rm = TRUE)
    }
  }
  
  c(Chi2 = sum(chi2, na.rm = TRUE))
}