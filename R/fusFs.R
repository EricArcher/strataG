#' @title Fu's Fs
#' @description Calculate Fu's Fs for a set of sequences to test 
#'   for selection.
#' 
#' @param x set of DNA sequences or a haploid \linkS4class{gtypes} 
#'   object with sequences.
#'   
#' @note Currently, this function is limited to calculating Fs for fewer than 
#'   approximately 172 sequences due to numerical overflow issues. \code{NaN} will 
#'   be returned for larger data sets.
#' 
#' @references Fu, Y-X. 1997. Statistical tests of neutrality of mutations 
#'   against population growth, hitchiking and background selection. 
#'   Genetics 147:915-925.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples 
#' data(dolph.seqs)
#' 
#' fusFs(dolph.seqs)
#' 
#' @importFrom copula Stirling1
#' @export
#' 
fusFs <- function(x) {
  fs.func <- function(h) {
    if(is.null(h)) return(NA)
    h$hap.seqs <- as.matrix(h$hap.seqs)
    h$haps <- na.omit(h$haps)
    n <- length(h$haps)
    k0 <- length(unique(h$haps))
    
    pws.diff <- dist.dna(h$hap.seqs, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
    pws.diff <- pws.diff[h$haps, h$haps]
    theta.pi <- mean(pws.diff[lower.tri(pws.diff)])
    # NaN produced for n > 172 from copula::Stirling1
    Sk.theta.k <- sapply(k0:n, function(k) abs(Stirling1(n, k)) * (theta.pi ^ k))
    # potential typo on Fu 1997 page 916 that implies prod(theta.pi - 0:(n + 1)). See Ewens 1972 Eqn 22.
    Sn.theta.pi <- prod(theta.pi + 0:(n - 1)) 
    s.prime <- sum(Sk.theta.k / Sn.theta.pi)
    log(s.prime / (1 - s.prime))
  }
  
  x <- if(inherits(x, "gtypes")) {
    sequences(x, as.haplotypes = FALSE)
  } else {
    as.multidna(x)
  }
  
  sapply(getSequences(x, simplify = FALSE), function(dna) {
    dna <- as.matrix(dna)
    rownames(dna) <- 1:nrow(dna)
    haps <- labelHaplotypes(dna)
    fs.func(haps)
  })
}
