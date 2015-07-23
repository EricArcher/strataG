#' @rdname popStructStat
#' @export
#' 
statFis <- function(g, strata = NULL, ...) {
  if(ploidy(g) < 2) return(c(Fis = NA))
  strata <- if(is.null(strata)) {
    strata(g)
  } else {
    rep(strata, length.out = nInd(g))
  }
  if(!is.factor(strata)) strata <- factor(strata)
  
  if(any(is.na(strata))) {
    toUse <- !is.na(strata)
    strata <- strata[toUse]
    g <- g[toUse, , ]
  }
  
  est <- statFis_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    as.numeric(strata) - 1,
    ploidy(g)
  )
  names(est) <- "Fis"
  est
}