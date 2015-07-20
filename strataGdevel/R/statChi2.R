#' @rdname popStructStat
#' @export
#' 
statChi2 <- function(g, strata = NULL, ...) {
  strata <- if(is.null(strata)) {
    strata(g)
  } else {
    rep(strata, length.out = nInd(g))
  }
  if(!is.factor(strata)) strata <- factor(strata)
  
  chi2 <- statChi2_C(
    sapply(loci(g), as.numeric),
    rep(as.numeric(factor(strata)), ploidy(g))
  )
  
#   chi2 <- vector("numeric", ncol(g@loci))
#   for(i in 1:length(chi2)) {
#     obs.freq <- table(g@loci[, i], strata, useNA = "no")
#     if(nrow(obs.freq) > 0 & ncol(obs.freq) > 1) {
#       exp.freq <- outer(rowSums(obs.freq), colSums(obs.freq)) / sum(obs.freq)
#       chi2[i] <- sum((obs.freq - exp.freq) ^ 2 / exp.freq, na.rm = TRUE)
#     }
#   }
  
  c(Chi2 = chi2)
}