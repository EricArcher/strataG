#' @rdname popStructStat
#' @export
#' 
statPhist <- function(g, strata = NULL, hap.dist = NULL, model = "K80", 
                      gamma = FALSE, pairwise.deletion = TRUE, ...)  {
  if(ploidy(g) != 1) return(c(PHIst = NA))
  stat.name <- "PHIst"
  
  # check hap.dist matrix
  if(!is.null(hap.dist)) {
    if(!("dist" %in% class(hap.dist))) {
      stop("'hap.dist' must be of class 'dist'")
    }
    if(all(hap.dist[lower.tri(hap.dist)] == 1)) stat.name <- "Fst"
  } else if(is.null(g@sequences)) {
    # format distances for Fst (all 1s, and 0s on the diagonal)
    haps <- unique(g@loci[, 1])
    hap.dist <- matrix(1, nrow = length(haps), ncol = length(haps), 
                       dimnames = list(haps, haps))
    diag(hap.dist) <- 0
    stat.name <- "Fst"
  } else {
    hap.dist <- dist.dna(
      getSequences(sequences(g), simplify = F)[[1]], model = model, gamma = gamma, 
      pairwise.deletion = pairwise.deletion
    )
  }
  if(!is.matrix(hap.dist)) hap.dist <- as.matrix(hap.dist)
  
  strata <- if(is.null(strata)) {
    strata(g)
  } else {
    rep(strata, length.out = nInd(g))
  }
  if(!is.factor(strata)) strata <- factor(strata)
  
  haps <- loci(g)[, 1]
  hap.dist <- hap.dist[levels(haps), levels(haps)]
  
  est <- statPhist_C(
    as.numeric(haps) - 1,
    as.numeric(strata) - 1,
    hap.dist
  )
  
#   # Extract summary values
#   strata.hap.freq <- table(g@loci[, 1], strata, useNA = "no")
#   haps <- rownames(strata.hap.freq)
#   hap.dist <- hap.dist[haps, haps, drop = FALSE]
#   strata.freq <- colSums(strata.hap.freq)
#   num.strata <- length(strata.freq)
#   num.samples <- sum(strata.freq)
#   
#   # Calculate sums of squares within strata (Eqn 8a)
#   ssd.wp <- sum(sapply(names(strata.freq), function(s) {
#     hap.freq <- strata.hap.freq[, s]
#     freq.prod <- outer(hap.freq, hap.freq)
#     sum(freq.prod * hap.dist) / (2 * strata.freq[s])
#   }))
#   
#   # Calculate sums of squares among strata (Eqn 8b)
#   pairs <- expand.grid(names(strata.freq), names(strata.freq))
#   ssd.ap <- sum(sapply(1:nrow(pairs), function(i) {
#     st <- unlist(pairs[i, ])
#     freq.prod <- outer(strata.hap.freq[, st[1]], strata.hap.freq[, st[2]])
#     sum(freq.prod * hap.dist)
#   })) / sum(2 * strata.freq)
#   ssd.ap <- ssd.ap - ssd.wp
#   
#   # Calculate average sample size correction for among strata variance 
#   #  Eqn 9a in paper, but modified as in Table 8.2.1.1 from Arlequin v3.5.1 manual
#   #  (denominator is sum{I} - 1)
#   n <- (num.samples - sum(strata.freq ^ 2 / num.samples)) / (num.strata - 1)
#   
#   # Calculate variance components (Table 1)
#   #   Set MSD (SSD / df) equal to expected MSD
#   Vc <- ssd.wp / (num.samples - num.strata)
#   Vb <- ((ssd.ap / (num.strata - 1)) - Vc) / n
#   
#   est <- Vb / (Vb + Vc)
  
  if(is.nan(est)) est <- NA
  names(est) <- stat.name
  est
}