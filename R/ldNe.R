#' @title ldNe
#' @description Estimate Ne from linkage disequilibrium based on Pearson correlation 
#'   approximation as in Waples et al 2016.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param maf.threshold smallest minimum allele frequency permitted to include 
#'   a locus in calculation.
#'
#' @return vector with sample size (\code{S}), number of pairwise loci comparisons 
#'   used (\code{num.comp}), mean r-squared (\code{mean.rsq}), and 
#'   Ne estimate (\code{Ne}).
#' 
#' @references Waples, R.S. 2006. A bias correction for estimates of effective population 
#'   size based on linkage disequilibrium at unlinked gene loci. 
#'   Conservation Genetics 7:167-184. \cr
#'   Waples RK, Larson WA, and Waples RS. 2016. 
#'   Estimating contemporary effective population size in non-model species using 
#'   linkage disequilibrium across thousands of loci. Heredity 117:233-240; 
#'   doi:10.1038/hdy.2016.60
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov} (adapted from code by R. Waples)
#' 
#' @export
#' 
ldNe <- function(g, maf.threshold = 0) {
  if(ploidy(g) != 2) stop("'g' must have diploid data")
  
  # remove non-biallelic loci
  num.alleles <- numAlleles(g)
  biallelic <- names(num.alleles)[num.alleles <= 2]
  if(length(biallelic) == 0) {
    warning("No loci are biallelic.")
    return(NULL)
  }
  g <- g[, biallelic, ]
  
  # remove loci below MAF threshold
  if(maf.threshold > 0) {
    maf.g <- maf(g)
    above.thresh <- which(maf.g >= maf.threshold)
    g <- g[, above.thresh, ]
  }
  
  # convert gtypes to coded data.frame
  df <- as.data.frame(g, ids = FALSE, strata = FALSE)
  new.g <- do.call(cbind, lapply(seq(1, ncol(df), by = 2), function(a1) {
    df.i <- df[, c(a1, a1 + 1)]
    alleles <- sort(unique(unlist(df.i)))
    apply(df.i, 1, function(i) sum(i == alleles[1]))
  }))
  
  # remove loci that are constant
  is.constant <- sapply(1:ncol(new.g), function(i) var(new.g[, i]) == 0)
  new.g <- new.g[, !is.constant]
  
  # calculate r-squared
  loc.pairs <- combn(ncol(new.g), 2)
  rsq <- sapply(1:ncol(loc.pairs), function(i) {
    l1 <- new.g[, loc.pairs[1, i]]
    l2 <- new.g[, loc.pairs[2, i]]
    cor(l1, l2, use = "pairwise.complete.obs", method = "pearson") 
  }) ^ 2
  
  # estimates Ne per LDNe
  S <- nrow(new.g)
  mean.rsq <- mean(rsq)
  Rprime <- mean.rsq * ((S / (S - 1)) ^ 2) - (1 / S) - (3.19 / S ^ 2)
  ne <- if(S > 29) {
    (1 / 3 + sqrt(1 / 9 - 2.76 * Rprime)) / (2 * Rprime)
  } else {
    (0.308 + sqrt(0.308 ^ 2 - 2.08 * Rprime)) / (2 * Rprime)
  }
  
  c(S = S, num.comp = ncol(loc.pairs), mean.rsq2 = mean.rsq, Ne = ne)
}