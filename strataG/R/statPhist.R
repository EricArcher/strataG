#' @rdname popStructStat
#' @importFrom apex getSequences
#' @importFrom swfscMisc harmonic.mean
#' @export
#' 
statPhist <- function(g, strata = NULL, hap.dist = NULL, model = "K80", 
                      gamma = FALSE, pairwise.deletion = TRUE, ...)  {
  if(ploidy(g) != 1 | nStrata(g) == 1) return(c(PHIst = NA))
  
  strata <- if(is.null(strata)) {
    strata(g)
  } else {
    rep(strata, length.out = nInd(g))
  }
  if(!is.factor(strata)) strata <- factor(strata)
  
  if(!is.null(hap.dist)) {
    if(inherits(hap.dist, "dist")) {
      hap.dist <- lapply(locNames(g), function(gene) hap.dist)
    } else if(is.list(hap.dist) & length(hap.dist) == nLoc(g)) {
      names(hap.dist) <- locNames(g)
    } else {
      warning("'hap.dist' is not NULL, nor a distance matrix or list as long as the number of genes.")
      hap.dist <- NULL
    }
  }
  
  est <- sapply(locNames(g), function(gene) {  
    stat.name <- "PHIst"
    # check hap.dist matrix
    hd <- hap.dist[[gene]]
    if(!is.null(hd)) {
      if(!("dist" %in% class(hd))) stop("'hap.dist' must be of class 'dist'")
      if(all(hd[lower.tri(hd)] == 1)) stat.name <- "Fst"
    } else if(is.null(sequences(g))) {
      # format distances for Fst (all 1s, and 0s on the diagonal)
      haps <- unique(as.character(loci(g)[, gene]))
      hd <- matrix(1, nrow = length(haps), ncol = length(haps), dimnames = list(haps, haps))
      diag(hd) <- 0
      stat.name <- "Fst"
    } else {
      hd <- dist.dna(
        getSequences(sequences(g, gene), simplify = TRUE), model = model, gamma = gamma, 
        pairwise.deletion = pairwise.deletion
      )
    }
    if(!is.matrix(hd)) hd <- as.matrix(hd)
    
    haps <- loci(g)[, gene]
    hd <- hd[levels(haps), levels(haps), drop = FALSE]
    
    result <- statPhist_C(as.numeric(haps) - 1, as.numeric(strata) - 1, hd)
    names(result) <- stat.name
    result
  }, USE.NAMES = FALSE)
  est[is.na(est)] <- 0
  stat.name <- if(any(names(est) == "PHIst")) "PHIst" else "Fst"
  est <- if(any(est <= 0)) mean(est, na.rm = TRUE) else harmonic.mean(est)
  
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