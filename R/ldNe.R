#' @title ldNe
#' @description Estimate Ne from linkage disequilibrium
#'
#' @param g a \linkS4class{gtypes} object.
#' @param maf.threshold 
#'
#' @return Ne estimate
#' 
#' @references Waples
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
ldNe <- function(g, maf.threshold = 0) {
  if(ploidy(g) != 2) stop("'g' must have diploid data")
  
  num.alleles <- numAlleles(g)
  biallelic <- names(num.alleles)[num.alleles <= 2]
  if(length(biallelic) == 0) {
    warning("No loci are biallelic.")
    return(NULL)
  }
  
  g <- g[, biallelic, ]
  new.g <- as.data.frame(g, one.col = TRUE, ids = FALSE, strata = FALSE)
  for(x in colnames(new.g)) new.g[, x] <- as.numeric(factor(new.g[, x])) - 1
  
  loc.pairs <- combn(ncol(new.g), 2)
  rsq <- sapply(1:ncol(loc.pairs), function(i) {
    l1 <- new.g[, loc.pairs[1, i]]
    l2 <- new.g[, loc.pairs[2, i]]
    # only calculate if both loci are not monomorphic
    if(var(l1, na.rm = TRUE) * var(l2, na.rm = TRUE) > 0) cor(l1, l2) else NA
  }) ^ 2
  
  # this accomplishes the same thing as the more complicated adjustment for sampling error in LDNe
  S <- nrow(new.g)
  Rprime <-  mean(rsq) - 1 / (S - 1)
  # estimates Ne per LDNe
  if(S > 29) {
    (2/3 + sqrt(4 / 9 - 7.2 * Rprime)) / (2 * Rprime)
  } else {
    (0.618 + sqrt(0.618 ^ 2 - 5.24 * Rprime)) / (2 * Rprime)
  }
}