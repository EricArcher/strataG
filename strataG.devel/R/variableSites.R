#' @title Variable Sites
#' @description Identify variable sites among sequences.
#'  
#' @param x a \code{\link{gtypes}} object with sequences, a list of sequences, or a consensus sequence.
#'   Sequences must be aligned. 
#' @param bases character vector of bases to consider.
#' 
#' @return A list with: \tabular{ll}{
#'   \code{site} \tab a list of sequences composed of variable sites. \cr
#'   \code{site.freqs} \tab a matrix of base pair frequencies by site. \cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{fixedSites}}
#' 
#' @export
#' 
variableSites <- function(x, bases = c("a", "c", "g", "t", "-")) {
  x <- as.multidna(x)
  result <- lapply(x@dna, function(dna) {
    site.freqs <- baseFreqs(dna, bases)$site.freqs
    var.site <- apply(site.freqs, 2, function(site) sum(site > 0) > 1)
    dna <- as.character(dna)
    var.seqs <- lapply(1:nrow(dna), function(i) dna[i, var.site])
    names(var.seqs) <- rownames(dna)
    var.site.freqs <- as.matrix(site.freqs[, var.site], nrow = nrow(site.freqs))
    colnames(var.site.freqs) <- which(var.site)
    list(sites = as.DNAbin(var.seqs), site.freqs = var.site.freqs)
  })
  
  if(length(result) == 1) {
    result[[1]]
  } else {
    names(result) <- names(x@dna)
    result
  }
}