#' @title Consensus Sequence
#' @description Return a consensus sequence from set of aligned sequences, 
#'   introducing IUPAC ambiguity codes where necessary.
#' 
#' @param x a \linkS4class{gtypes} object with aligned sequences or a list of 
#'   aligned DNA sequences.
#' @param ignore.gaps logical. Ignore gaps at a site when creating consensus. 
#'   If \code{TRUE}, then bases with a gap are removed before consensus is 
#'   calculated. If \code{FALSE} and a gap is present, then the result is a gap.
#' @param simplify if there is a single locus, return result in a simplified
#'   form? If \code{FALSE} a list will be returned wth one element per locus.
#'   
#' @return A character vector of the consensus sequence.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.seqs)
#' createConsensus(dolph.seqs)
#' 
#' @export
#' 
createConsensus <- function(x, ignore.gaps = FALSE, simplify = TRUE) { 
  result <- sapply(
    apex::getSequences(as.multidna(x, as.haplotypes = FALSE), simplify = FALSE), 
    function(dna) {
      dna <- as.character(as.matrix(dna))
      apply(dna, 2, iupacCode, ignore.gaps = ignore.gaps)
    },
    simplify = FALSE
  )
  
  if(length(result) == 1 & simplify) result[[1]] else result
}