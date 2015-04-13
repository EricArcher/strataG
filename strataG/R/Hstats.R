#' @title Heterozygosity estimates: Ho, Hs, and Ht
#'
#' @description Calculate overall heterozygosity (Ho), expected heterozygosity in
#'   subpopulations (Hs), and expected total heterozygosity (Ht)
#'
#' @param g a \code{gtypes} object
#'
#' @return matrix of Ho, Hs, and Ht (rows) for each locus (columns)
#'
#' @author Eric Archer <eric.archer@@noaa.gov>
#'
#' @references Nei, M. and R.K. Chesser. 1983. Estimation of fixation indices and gene diversities. Ann. Hum. Genet. 47:253-259.
#'
#' @export

Hstats <- function(g) {
  #stopifnot.gtypes(g, "diploid")

  nloc <- ncol(g@loci)
  result <- matrix(0, nrow = 3, ncol = nloc)
  rownames(result) <- c("Ho", "Hs", "Ht")
  colnames(result) <- colnames(g@loci)
  strata <- rep(g@strata, g@ploidy)
  for(i in 1:nloc) {
    ## Estimate Ho (frequency of all heterozygotes): Equation 5, page 254
    locus <- g@loci[, i]
    loc.mat <- matrix(locus, ncol = g@ploidy)
    genotype <- sapply(1:nrow(loc.mat), function(i) {
      alleles <- loc.mat[i, ]
      if(any(is.na(alleles))) return(NA)
      if(length(unique(alleles)) == 1) alleles[1] else "het"
    })
    hom.freq <- prop.table(table(g@strata, genotype), 1)
    hom.freq <- hom.freq[, colnames(hom.freq) != "het", drop = FALSE]
    Ho <- if(ncol(hom.freq) == 0) 0 else {
      mean.hom.freq <- colMeans(hom.freq)
      1 - sum(mean.hom.freq)
    }

    ## Estimate Hs (expected heterozygosity within strata): Equation 9, page 255
    allele.freq <- table(strata, locus)
    strata.freq <- rowSums(allele.freq) / g@ploidy
    hom.freq <- prop.table(allele.freq, 1)
    mean.het <- mean(1 - rowSums(hom.freq ^ 2))
    harm.n <- harmonic.mean(strata.freq)
    Hs <- (harm.n / (harm.n - 1)) * (mean.het - (Ho / 2 / harm.n))

    ## Estimate Ht (expected heterozygosity overall): Equation 11, page 256
    mean.het <- 1 - sum(colMeans(hom.freq) ^ 2)
    harm.n.s <- harm.n * sum(allele.freq)
    Ht <- mean.het + (Hs / harm.n.s) - (Ho / 2 / harm.n.s)
    result[, i] <- c(Ho, Hs, Ht)
  }
  result
}