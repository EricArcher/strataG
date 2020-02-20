#' @title ldNe
#' @description Estimate Ne from linkage disequilibrium based on Pearson 
#'   correlation approximation following Waples et al 2016. Adapted from code 
#'   by R. Waples and W. Larson.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param maf.threshold smallest minimum allele frequency permitted to include 
#'   a locus in calculation of Ne.
#' @param by.strata apply the `maf.threshold` by strata. If `TRUE`
#'   then loci that are below this threshold in any strata will be removed 
#'   from the calculation of Ne for all strata. Loci below `maf.threshold` 
#'   within a stratum are always removed for calculations of Ne for that 
#'   stratum.
#' @param ci central confidence interval.
#' @param drop.missing drop loci with missing genotypes? If `FALSE`, a slower 
#'   procedure is used where individuals with missing genotypes are removed 
#'   in a pairwise fashion. 
#' @param num.cores The number of cores to use to distribute computations over.
#'   If set to \code{NULL}, the value will be what is reported 
#'   by \code{\link[parallel]{detectCores} - 1}.
#'
#' @return a data.frame with one row per strata and the following columns:
#' \describe{
#'  \item{\code{stratum}}{stratum being summarized}
#'  \item{\code{S}}{harmonic mean of sample size across pairwise comparisons of
#'  loci}
#'  \item{\code{num.comp}}{number of pairwise loci comparisons used}
#'  \item{\code{mean.rsq}}{mean r^2 over all loci}
#'  \item{\code{mean.E.rsq}}{mean expected r^2 over all loci}
#'  \item{\code{Ne}}{estimated Ne}
#'  \item{\code{param.lci, param.uci}}{parametric lower and upper CIs}
#' }
#' 
#' @references Waples, R.S. 2006. A bias correction for estimates of effective
#'   population size based on linkage disequilibrium at unlinked gene loci.
#'   Conservation Genetics 7:167-184. \cr
#'   Waples RK, Larson WA, and Waples RS. 2016. Estimating contemporary 
#'   effective population size in non-model species using linkage 
#'   disequilibrium across thousands of loci. Heredity 117:233-240; 
#'   doi:10.1038/hdy.2016.60
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
ldNe <- function(g, maf.threshold = 0, by.strata = FALSE, ci = 0.95, 
                 drop.missing = FALSE, num.cores = 1) {
  if(getPloidy(g) != 2) stop("'g' must have diploid data")
  
  mat <- as.data.frame(g, coded.snps = TRUE) %>% 
    tibble::column_to_rownames("id") 
  st <- mat$stratum
  mat$stratum <- NULL
  mat <- as.matrix(mat)
  
  # remove loci if below maf threshold for any stratum
  if(maf.threshold > 0 & by.strata) {
    above.thresh <- do.call(cbind, tapply(1:nrow(mat), st, function(i) {
      colMeans(mat[i, ]) / 2
    })) %>% 
      apply(1, function(x) all(x >= maf.threshold)) %>% 
      which() %>% 
      names()
    if(length(above.thresh) < 2) {
      warning(
        "Fewer than two loci are above 'maf.threshold' in all strata.",
        "NULL returned.", call. = FALSE
      )
      return(NULL)
    }
    mat <- mat[, above.thresh, ]
  }
  
  # calculate Pearson r-squared between a pair of loci
  .compLoc <- function(i, loc.pairs, mat) {
    pair.mat <- mat[, loc.pairs[, i]]
    pair.mat <- pair.mat[stats::complete.cases(pair.mat), , drop = FALSE]
    rsq <- stats::cor(pair.mat[, 1], pair.mat[, 2], method = "pearson") ^ 2
    S <- nrow(pair.mat)
    c(S = S, rsq = rsq)
  }
  
  .calcRsqMissing <- function(mat) {
    # calculate correlation r-squared (rsq) between all pairs of loci
    loc.pairs <- utils::combn(ncol(mat), 2)
    cl <- swfscMisc::setupClusters(num.cores)
    loc.comp.mat <- tryCatch({
      if(is.null(cl)) {
        lapply(1:ncol(loc.pairs), .compLoc, loc.pairs = loc.pairs, mat = mat)
      } else {
        parallel::clusterEvalQ(cl, require(strataG))
        parallel::clusterExport(cl, c("loc.pairs", "mat"), environment())
        parallel::parLapply(
          cl, 1:ncol(loc.pairs), .compLoc, loc.pairs = loc.pairs, mat = mat
        )
      } 
    }, finally = if(!is.null(cl)) parallel::stopCluster(cl))
    do.call(cbind, loc.comp.mat)
    # to.keep <- apply(loc.comp.mat, 2, function(x) all(!is.na(x)))
    # if(sum(to.keep) == 0) stop("all loci have missing data")
    # loc.comp.mat[, to.keep, drop = FALSE]
  }
  
  .calcRsq <- function(x) {
    opts <- options(matprod = if(any(is.na(x))) "default" else "blas")
    # Transpose to make matrix math work
    x <- t(x)
    # Center each variable
    x <- x - rowMeans(x, na.rm = TRUE)
    # Standardize each variable
    x <- x / sqrt(rowSums(x ^ 2, na.rm = TRUE))
    # Calculate correlations
    rsq <- tcrossprod(x) ^ 2 
    options(opts)
    mean(rsq[lower.tri(rsq)])
  }
  
  # Eqn 1.7: calculate Ne
  .calcNe <- function(S, Rsq.drift) {
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
  ne.smry <- tapply(1:nrow(mat), st, function(i) {
    mat.st <- mat[i, ]
    
    # remove loci below MAF threshold
    if(maf.threshold > 0) {
      above.thresh <- (colMeans(mat.st) / 2) >= maf.threshold
      which() %>% 
        names()
      if(length(above.thresh) < 2) {
        warning(
          "Fewer than two loci are above 'maf.threshold' in", 
          paste0("'", unique(st[i]), "'."),
          "NULL returned.", call. = FALSE
        )
        return(NULL)
      }
      mat.st <- mat.st[, above.thresh]
    }
    
    # keep only polymorphic loci
    polymorph <- apply(mat.st, 2, function(x) length(unique(x)) > 1)
    if(sum(polymorph) < 2) {
      warning(
        "Fewer than two loci are polymorphic in", 
        " '", unique(st[i]), "'. NULL returned.", 
        call. = FALSE
      )
      return(NULL)
    }
    mat.st <- mat.st[, polymorph, drop = FALSE]
    
    # calculate r-squared among pairs of loci
    # use matrix algebra if no missing data
    loc.missing <- which(apply(mat.st, 2, function(x) any(is.na(x))))
    rsq.list <-  if(length(loc.missing) == 0) {
      rsq <- .calcRsq(mat.st)
      N <- (ncol(mat.st) * (ncol(mat.st) - 1)) / 2
      list(rsq = rsq, S = nrow(mat.st), N = N)
    } else if(drop.missing) {
      mat.st <- mat.st[, -loc.missing, drop = FALSE]
      if(ncol(mat.st) >= 2) {
        rsq <- .calcRsqMissing(mat.st)
        list(rsq = rsq["rsq", ], S = rsq["S", ], N = ncol(rsq))
      } else {
        warning(
          "Can't compute ldNe in '", unique(st[i]), "' ",
          "because fewer than 2 loci are missing genotypes. NULL returned.", 
          call. = FALSE
        )
        return(NULL)
      }
    } else {
      warning(
        "Can't compute ldNe in '", unique(st[i]), "' ",
        "because loci are missing genotypes and 'drop.missing = FALSE'. ",
        "NULL returned.", call. = FALSE
      )
      return(NULL)
    }
    
    S <- rsq.list$S
    # Eqn 1.1: expected r-squared
    E.rsq <- ifelse(
      S > 29, 
      (1 / S) + (3.19 / S ^ 2), 
      0.0018 + (0.907 / S) + (4.44 / S ^ 2)
    )
    # sample size corrected r-squared
    rsq <- rsq.list$rsq * ((S / (S - 1)) ^ 2)
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
    N <- rsq.list$N
    S.harm.mean <- N / sum(1 / S)
    # initial Ne.0
    ne0 <- .calcNe(S.harm.mean, Rsq.drift.0)
    
    # Eqn 1.11: R-squared prime weights
    wt <- S ^ 2 / (S + 3 * ne0) ^ 2
    # Eqn 1.12: weighted R-squared prime
    Rsq.drift <- sum(wt * Rsq.drift) / sum(wt)
    
    # mean expected r-squared
    mean.E.rsq <- sum(w * E.rsq) / W
    
    # calculate CI
    lci.p <- (1 - ci) / 2
    uci.p <- 1 - lci.p
    Rsq.drift.lci <- mean.rsq * N / stats::qchisq(lci.p, N) - mean.E.rsq
    Rsq.drift.uci <- mean.rsq * N / stats::qchisq(uci.p, N) - mean.E.rsq
    
    ne <- .calcNe(S.harm.mean, Rsq.drift)
    param.lci <- .calcNe(S.harm.mean, Rsq.drift.lci)
    param.uci <- .calcNe(S.harm.mean, Rsq.drift.uci)
    if(ne < 0) ne <- Inf
    if(param.lci < 0) param.lci <- Inf
    if(param.uci < 0) param.uci <- Inf
    
    c(
      S = S.harm.mean, num.comp = N, mean.rsq = mean.rsq, 
      mean.E.rsq = mean.E.rsq, Ne = ne, param.lci = param.lci, 
      param.uci = param.uci
    )
  })
  
  ne.smry <- ne.smry[!sapply(ne.smry, is.null)]
  if(length(ne.smry) == 0) return(NULL)
  do.call(rbind, ne.smry) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("stratum") %>% 
    dplyr::select(.data$stratum, dplyr::everything())
}