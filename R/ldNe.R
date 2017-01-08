#' @title ldNe
#' @description Estimate Ne from linkage disequilibrium based on Pearson 
#'   correlation approximation following Waples et al 2016. Adapted from code 
#'   by R. Waples and W. Larson.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param maf.threshold smallest minimum allele frequency permitted to include 
#'   a locus in calculation of Ne.
#' @param by.strata apply the \code{maf.threshold} by strata. If \code{TRUE}
#'   then any locus that is below this threshold in any strata will be removed 
#'   from the calculation of Ne for every stratum. Otherwise, loci are removed 
#'   only if they are below the \code{maf.threshold} in the stratum for which 
#'   Ne is calculated.
#' @param ci central confidence interval.
#'
#' @return a numeric matrix with one row per strata and the following columns:
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
#'   Waples RK, Larson WA, and Waples RS. 2016. Estimating contemporary 
#'   effective population size in non-model species using linkage 
#'   disequilibrium across thousands of loci. Heredity 117:233-240; 
#'   doi:10.1038/hdy.2016.60
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom parallel mclapply detectCores
#' @importFrom stats cor qchisq
#' @export
#' 
ldNe <- function(g, maf.threshold = 0, by.strata = FALSE, ci = 0.95) {
  if(ploidy(g) != 2) stop("'g' must have diploid data")
  
  # remove non-biallelic loci
  num.alleles <- numAlleles(g)
  biallelic <- names(num.alleles)[num.alleles <= 2]
  if(length(biallelic) == 0) {
    warning("No loci are biallelic. NULL returned.", call. = FALSE)
    return(NULL)
  }
  g <- g[, biallelic, ]
  
  # remove loci if below maf threshold for any stratum
  if(maf.threshold > 0 & by.strata) {
    maf.g <- maf(g, by.strata = TRUE)
    above.thresh <- which(apply(maf.g, 2, function(x) all(x >= maf.threshold)))
    if(length(above.thresh) < 2) {
      warning(
        paste0("Fewer than two loci are above 'maf.threshold' in any stratum. NULL returned."),
        call. = FALSE
      )
      return(NULL)
    }
    g <- g[, above.thresh, ]
  }
  
  # calculate Pearson r-squared between a pair of loci
  compLoc <- function(i, loc.pairs, mat) {
    pair.mat <- mat[, loc.pairs[, i]]
    pair.mat <- pair.mat[complete.cases(pair.mat), , drop = FALSE]
    rsq <- cor(pair.mat[, 1], pair.mat[, 2], method = "pearson") ^ 2
    S <- nrow(pair.mat)
    c(S = S, rsq = rsq)
  }
  
  # Eqn 1.7: calculate Ne
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
  
  # calculate Ne by strata
  ne.smry <- sapply(strataSplit(g), function(g.st) {
    # remove loci below MAF threshold
    if(maf.threshold > 0) {
      above.thresh <- which(maf(g.st) >= maf.threshold)
      if(length(above.thresh) < 2) {
        warning(
          paste0("Fewer than two loci are above 'maf.threshold' in '", strataNames(g.st), "'"),
          call. = FALSE
        )
        return(NULL)
      }
      g.st <- g.st[, above.thresh, ]
    }
    
    # create list of coded numeric matrices for each stratum
    g.mat <- as.matrix(g.st, ids = FALSE, strata = FALSE)
    mat <- matrix(nrow = nInd(g.st), ncol = nLoc(g.st))
    a1.i <- seq(1, ncol(g.mat), by = 2)
    for(j in 1:length(a1.i)) {
      a1 <- a1.i[j]
      mat.j <- g.mat[, c(a1, a1 + 1), drop = FALSE]
      alleles <- sort(unique(mat.j))
      for(i in 1:nrow(mat.j)) mat[i, j] <- sum(mat.j[i, ] == alleles[1])
    }
    
    # remove loci that are constant
    mat <- mat[, apply(mat, 2, function(x) var(x) > 0), drop = FALSE]
    if(ncol(mat) < 2) {
      warning(
        paste0("Fewer than two loci have more than one genotype in '", strataNames(g.st), "'"),
        call. = FALSE
      )
      return(NULL)
    }
    
    # calculate correlation r-squared (rsq) between pairs of loci
    loc.pairs <- combn(ncol(mat), 2)
    loc.comp.mat <- mclapply(1:ncol(loc.pairs), compLoc, loc.pairs = loc.pairs, mat = mat, mc.cores = detectCores() - 1)
    loc.comp.mat <- do.call(cbind, loc.comp.mat)
    # cl <- .setupClusters()
    # loc.comp.mat <- tryCatch({
    #   if(!is.null(cl)) {
    #     parSapply(cl, 1:ncol(loc.pairs), compLoc, loc.pairs = loc.pairs, mat = mat)
    #   } else {
    #     sapply(1:ncol(loc.pairs), compLoc, loc.pairs = loc.pairs, mat = mat)
    #   }
    # }, finally = if(!is.null(cl)) stopCluster(cl))
    
    S <- loc.comp.mat["S", ]
    # Eqn 1.1: expected r-squared
    E.rsq <- ifelse(S > 29, (1 / S) + (3.19 / S ^ 2), 0.0018 + (0.907 / S) + (4.44 / S ^ 2))
    # sample size corrected r-squared
    rsq <- loc.comp.mat["rsq", ] * ((S / (S - 1)) ^ 2)
    # Eqn 1.4
    w <- S ^ 2
    # Eqn 1.5
    W <- sum(w)
    mean.rsq <- sum(w * rsq) / W
    # Eqn 1.6
    Rsq.drift <- rsq - E.rsq
    # Eqn 1.10: R-squared prime.0 for Ne.0 
    Rsq.drift.0 <- sum(Rsq.drift * w) / W
    # harmonic mean of S
    N <- ncol(loc.comp.mat)
    S.harm.mean <- N / sum(1 / S)
    # initial Ne.0
    ne0 <- calcNe(S.harm.mean, Rsq.drift.0)
    
    # Eqn 1.11: R-squared prime weights
    wt <- S ^ 2 / (S + 3 * ne0) ^ 2
    # Eqn 1.12: weighted R-squared prime
    Rsq.drift <- sum(wt * Rsq.drift) / sum(wt)
    
    # mean expected r-squared
    mean.E.rsq <- sum(w * E.rsq) / W
    
    # calculate CI
    lci.p <- (1 - ci) / 2
    uci.p <- 1 - lci.p
    Rsq.drift.lci <- mean.rsq * N / qchisq(lci.p, N) - mean.E.rsq
    Rsq.drift.uci <- mean.rsq * N / qchisq(uci.p, N) - mean.E.rsq
    
    ne <- calcNe(S.harm.mean, Rsq.drift)
    param.lci <- calcNe(S.harm.mean, Rsq.drift.lci)
    param.uci <- calcNe(S.harm.mean, Rsq.drift.uci)
    if(ne < 0) ne <- Inf
    if(param.lci < 0) param.lci <- Inf
    if(param.uci < 0) param.uci <- Inf
    
    c(
      S = S.harm.mean, num.comp = N, mean.rsq = mean.rsq, 
      mean.E.rsq = mean.E.rsq, Ne = ne, param.lci = param.lci, 
      param.uci = param.uci
    )
  }, USE.NAMES = TRUE, simplify = FALSE)
  ne.smry <- ne.smry[!sapply(ne.smry, is.null)]
  do.call(rbind, ne.smry)
}