#' @title Variable Sites
#' @description Identify variable sites among sequences.
#'  
#' @param x a \linkS4class{gtypes} object with sequences, 
#'   a \code{\link[ape]{DNAbin}} object, or a list of sequences.
#' @param bases character vector of bases to consider.
#' @param simplify if there is a single locus, return result in a simplified
#'   form? If \code{FALSE} a list will be returned wth one element per locus.
#' 
#' @return A list with: \describe{
#'   \item{site}{a \code{\link[ape]{DNAbin}} object composed of variable sites.}
#'   \item{site.freqs}{a matrix of base pair frequencies by site.}
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.haps)
#' 
#' variableSites(dolph.haps)
#' 
#' @seealso \code{\link{fixedSites}}
#' 
#' @export
#' 
variableSites <- function(x, bases = c("a", "c", "g", "t", "-"), simplify = TRUE) {
  bases <- tolower(bases)
  result <- sapply(
    apex::getSequences(as.multidna(x), simplify = FALSE),
    function(dna) {
      dna <- as.matrix(dna)
      site.freqs <- baseFreqs(dna, bases)$site.freqs
      var.site <- apply(site.freqs, 2, function(site) sum(site > 0) > 1)
      dna <- as.character(dna)
      var.seqs <- stats::setNames(
        lapply(1:nrow(dna), function(i) dna[i, var.site]),
        rownames(dna)
      )
      site.freqs <- site.freqs[, var.site, drop = FALSE]
      colnames(site.freqs) <- which(var.site)
      list(
        sites = as.matrix(ape::as.DNAbin(var.seqs)), 
        site.freqs = site.freqs
      )
    },
    simplify = FALSE
  )
  
  if(length(result) == 1 & simplify) result[[1]] else result
}