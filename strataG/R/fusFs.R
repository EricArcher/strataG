#' @title Fu's Fs
#' @description Calculate Fu's Fs for a set of sequences to test 
#'   for selection.
#' 
#' @param x set of DNA sequences or a haploid \linkS4class{gtypes} 
#'   object with sequences.
#' 
#' @references Fu, Y-X. 1997. Statistical tests of neutrality of mutations 
#'   against population growth, hitchiking and background selection. 
#'   Genetics 147:915-925.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom copula Stirling1
#' @importFrom apex getSequences
#' @export
#' 
fusFs <- function(x) {
  x <- as.multidna(x)
  
  sapply(getSequences(x, simplify = FALSE), function(dna) {
    dna <- as.matrix(dna)
    haps <- as.matrix(labelHaplotypes(dna)$hap.seqs)
    n <- nrow(dna)
    k0 <- nrow(haps)
    
    pws.diff <- dist.dna(haps, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
    theta.pi <- mean(pws.diff[lower.tri(pws.diff)])
    Sn.theta.pi <- prod(theta.pi - 0:(n + 1))
    Sk.theta.k <- sapply(k0:n, function(k) abs(Stirling1(n, k)) * (theta.pi ^ k))
    s.prime <- sum(Sk.theta.k / Sn.theta.pi)
    
    log(s.prime / (1 - s.prime))
  })
}
