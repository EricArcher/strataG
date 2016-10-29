#' @title ldNe
#' @description Estimate Ne from linkage disequilibrium based on Pearson 
#'   correlation approximation following Waples et al 2016. Adapted from code 
#'   by R. Waples and W. Larson.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param maf.threshold smallest minimum allele frequency permitted to include 
#'   a locus in calculation of Ne.
#' @param by.strata apply the \code{maf.threshold} by strata. if \code{TRUE} 
#'   then any locus that is below this threshold in any strata will be removed 
#'   from the calculation of Ne.
#' @param ci central confidence interval.
#'
#' @return a named numeric vector with: 
#' \describe{
#'  \item{\code{S}}{harmonic mean of sample size across pairwise comparisons of loci}
#'  \item{\code{num.comp}}{number of pairwise loci comparisons used}
#'  \item{\code{mean.rsq}}{mean r^2 over all loci}
#'  \item{\code{mean.E.rsq}}{mean expected r^2 over all loci}
#'  \item{\code{Ne}}{estimated Ne}
#'  \item{\code{param.lci, param.uci}}{parametric lower and upper CIs}
#' }
#' 
#' @references Waples, R.S. 2006. A bias correction for estimates of effective population 
#'   size based on linkage disequilibrium at unlinked gene loci. 
#'   Conservation Genetics 7:167-184. \cr\cr
#'   Waples RK, Larson WA, and Waples RS. 2016. 
#'   Estimating contemporary effective population size in non-model species using 
#'   linkage disequilibrium across thousands of loci. Heredity 117:233-240; 
#'   doi:10.1038/hdy.2016.60
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom parallel mclapply
#' @importFrom stats cor qchisq
#' @export
#' 
ldNe <- function(g, maf.threshold = 0, by.strata = FALSE, ci = 0.95) {
  if(ploidy(g) != 2) stop("'g' must have diploid data")
  
  # remove non-biallelic loci
  num.alleles <- numAlleles(g)
  biallelic <- names(num.alleles)[num.alleles <= 2]
  if(length(biallelic) == 0) {
    warning("No loci are biallelic. NULL returned.")
    return(NULL)
  }
  g <- g[, biallelic, ]
  
  # remove loci below MAF threshold
  if(maf.threshold > 0) {
    maf.g <- maf(g, by.strata = by.strata)
    above.thresh <- if(by.strata) {
      which(apply(maf.g, 2, function(x) all(x >= maf.threshold)))
    } else {
      which(maf.g >= maf.threshold)
    }
    g <- g[, above.thresh, ]
  }
  
  # calculate Ne by strata
  t(sapply(strataSplit(g), function(st.g) {
    # convert gtypes to coded data.frame
    g.df <- as.data.frame(st.g, ids = FALSE, strata = FALSE)
    df <- do.call(cbind, lapply(seq(1, ncol(g.df), by = 2), function(a1) {
      df.i <- g.df[, c(a1, a1 + 1)]
      alleles <- sort(unique(unlist(df.i)))
      apply(df.i, 1, function(i) sum(i == alleles[1]))
    }))
    
    # remove loci that are constant
    is.constant <- sapply(1:ncol(df), function(i) var(df[, i]) == 0)
    df <- df[, !is.constant]
    
    # Eqn 1.7
    calcNe <- function(S, Rsq.drift) {
      if(S > 29) {
        root.term <- 1 / 9 - 2.76 * Rsq.drift
        if(root.term < 0) root.term <- 0
        (1 / 3 + sqrt(root.term)) / (2 * Rsq.drift)
      } else {
        root.term <- 0.308 ^ 2 - 2.08 * Rsq.drift
        if(root.term < 0) root.term <- 0
        (0.308 + sqrt(root.term)) / (2 * Rsq.drift)
      }
    }
    
    # calculate correlation r-squared (rsq) between pairs of loci
    loc.pairs <- combn(ncol(df), 2)
    num.cores <- if(.Platform$OS.type == "windows") 1 else detectCores() - 1
    loc.comp.mat <- do.call(rbind, mclapply(1:ncol(loc.pairs), function(i) {
      pair.df <- df[, loc.pairs[, i]]
      pair.df <- pair.df[complete.cases(pair.df), ]
      rsq <- cor(pair.df[, 1], pair.df[, 2], method = "pearson") ^ 2
      S <- nrow(pair.df)
      c(S = S, rsq = rsq)
    }, mc.cores = num.cores ))
    
    S <- loc.comp.mat[, "S"]
    # Eqn 1.1: expected r-squared
    E.rsq <- ifelse(S > 29, (1 / S) + (3.19 / S ^ 2), 0.0018 + (0.907 / S) + (4.44 / S ^ 2))
    # sample size corrected r-squared
    rsq <- loc.comp.mat[, "rsq"] * ((S / (S - 1)) ^ 2)
    # Eqn 1.4
    w <- loc.comp.mat[, "S"] ^ 2
    # Eqn 1.5
    W <- sum(w)
    mean.rsq <- sum(w * rsq) / W
    # Eqn 1.6
    Rsq.drift <- rsq - E.rsq
    # Eqn 1.10: R-squared prime.0 for Ne.0 
    Rsq.drift.0 <- sum(Rsq.drift * w) / W
    # harmonic mean of S
    S.harm.mean <- nrow(loc.comp.mat) / sum(1 / S)
    # initial Ne.0
    ne0 <- calcNe(S.harm.mean, Rsq.drift.0)
    
    # Eqn 1.11: R-squared prime weights
    wt <- S ^ 2 / (S + 3 * ne0) ^ 2
    # Eqn 1.12: weighted R-squared prime
    Rsq.drift <- sum(wt * Rsq.drift) / sum(wt)
    
    # mean expected r-squared
    mean.E.rsq <- sum(w * E.rsq) / W
    
    # calculate CI
    N <- nrow(loc.comp.mat)
    lci.p <- (1 - ci) / 2
    uci.p <- 1 - lci.p
    Rsq.drift.lci <- mean.rsq * N / qchisq(lci.p, N) - mean.E.rsq
    Rsq.drift.uci <- mean.rsq * N / qchisq(uci.p, N) - mean.E.rsq
    
    c(
      S = S.harm.mean, 
      num.comp = N, 
      mean.rsq = mean.rsq,
      mean.E.rsq = mean.E.rsq,
      Ne = calcNe(S.harm.mean, Rsq.drift),
      param.lci = calcNe(S.harm.mean, Rsq.drift.lci),
      param.uci = calcNe(S.harm.mean, Rsq.drift.uci)
    )
  }))
}