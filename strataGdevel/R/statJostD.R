#' @rdname popStructStat
#' @export
#' 
statJostD <- function(g, strata = NULL, ...) {
  if(ploidy(g) == 1) return(c(D = NA))
  
  strata <- if(is.null(strata)) {
    strata(g)
  } else {
    rep(strata, length.out = nInd(g))
  }
  if(!is.factor(strata)) strata <- factor(strata)
  
  est <- statJostD_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    as.numeric(strata) - 1,
    ploidy(g)
  )
  
#   terms <- sapply(1:ncol(g@loci), function(i) {
#     allele.strata.freq <- table(g@loci[, i], strata, useNA = "no")
#     num.strata <- ncol(allele.strata.freq)
#     i.terms <- sapply(rownames(allele.strata.freq), function(allele) {
#       j.terms <- sapply(colnames(allele.strata.freq), function(strata) {
#         Nj <- sum(allele.strata.freq[, strata])
#         Nij <- allele.strata.freq[allele, strata]
#         a.term1 <- Nij / Nj
#         a.term2 <- a.term1 ^ 2
#         b.term <- Nij * (Nij - 1) / (Nj * (Nj - 1))
#         c(a.term1 = a.term1, a.term2 = a.term2, b.term = b.term)
#       })
#       a.term1 <- sum(j.terms["a.term1", ]) ^ 2
#       a.term2 <- sum(j.terms["a.term2", ])
#       c(a = (a.term1 - a.term2) / (num.strata - 1), 
#         b = sum(j.terms["b.term", ])
#       )
#     })
#     
#     c(a = sum(i.terms["a", ]), b = sum(i.terms["b", ]))
#   })
#   est <- harmonic.mean(1 - terms["a", ] / terms["b", ])

  c(D = est)
}