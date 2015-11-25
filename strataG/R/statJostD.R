#' @rdname popStructStat
#' @importFrom swfscMisc harmonic.mean 
#' @export
#' 
statJostD <- function(g, strata = NULL, ...) {
  if(ploidy(g) < 2 | nStrata(g) == 1) return(c(D = NA))
  
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
  
#   est <- statJostD_C(
#     sapply(loci(g), function(x) as.numeric(x) - 1), 
#     as.numeric(strata) - 1,
#     ploidy(g)
#  )
  
  allele.freqs <- alleleFreqs(g, by.strata = TRUE)
  terms <- sapply(locNames(g), function(x) {
    allele.strata.freq <- allele.freqs[[x]][, "freq", ]
    num.strata <- ncol(allele.strata.freq)
    i.terms <- sapply(rownames(allele.strata.freq), function(a) {
      j.terms <- sapply(colnames(allele.strata.freq), function(s) {
        Nj <- sum(allele.strata.freq[, s])
        Nij <- allele.strata.freq[a, s]
        a.term1 <- Nij / Nj
        a.term2 <- a.term1 ^ 2
        b.term <- Nij * (Nij - 1) / (Nj * (Nj - 1))
        c(a.term1 = a.term1, a.term2 = a.term2, b.term = b.term)
      })
      a.term1 <- sum(j.terms["a.term1", ]) ^ 2
      a.term2 <- sum(j.terms["a.term2", ])
      c(a = (a.term1 - a.term2) / (num.strata - 1), 
        b = sum(j.terms["b.term", ])
      )
    })
    
    c(a = sum(i.terms["a", ]), b = sum(i.terms["b", ]))
  })
  d.by.locus <- 1 - terms["a", ] / terms["b", ]
  d.by.locus <- ifelse(d.by.locus < 0, 0, d.by.locus)
  est <- harmonic.mean(d.by.locus)

  c(D = est)
}