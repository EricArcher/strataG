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
#' @export
#' 
FusFs <- function(x) {
  x <- as.multidna(x)
  
  result <- do.call(rbind, lapply(getSequences(x, simplify = FALSE), function(dna) {
    dna <- as.matrix(labelHaplotypes(dna)$hap.seqs)
    num.sites <- ncol(dna)
    
    pws.diff <- dist.dna(dna, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
    theta.pi <- mean(pws.diff[lower.tri(pws.diff)])
    
    Sn.theta.n <- prod(theta.pi - 0:(ncol(dna) + 1))
    
    s.prime <- 1
    
    log(s.prime / (1 - s.prime))
  }))
  
  if(nrow(result) == 1) {
    result[1, ]
  } else {
    rownames(result) <- locusNames(x)
    result
  }
}
