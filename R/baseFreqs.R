#' @title Base Frequencies
#' @description Calculate nucleotide base frequencies along a sequence.
#'
#' @param x a \linkS4class{gtypes} object with aligned sequences or a list of
#'   aligned DNA sequences.
#' @param bases character vector of bases. Must contain valid IUPAC codes. If
#'   \code{NULL}, will return summary of frequencies of observed bases.
#' @param ignore a character vector of bases to ignore when calculating site
#'   frequencies.
#' @param simplify if there is a single locus, return result in a simplified
#'   form? If \code{FALSE} a list will be returned wth one element per locus.
#'
#' @return For each gene, a list containing: 
#' \tabular{ll}{ 
#'   \code{site.freqs} \tab a matrix of base frequencies at each site.\cr 
#'   \code{base.freqs} \tab a vector of overall base proportion composition.\cr 
#' }
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(dloop.g)
#' bf <- baseFreqs(dloop.g)
#'
#' # Frequencies of first 10 sites
#' bf$site.freqs[, 1:10]
#'
#' # Base composition
#' bf$base.freqs
#'
#' @export
#' 
baseFreqs <- function(x, bases = NULL, ignore = c("n", "x", "-", "."), 
                      simplify = TRUE) {
  bases <- if(is.null(bases)) {
    rownames(iupac.mat)
  } else {
    tolower(as.character(bases))
  }
  ignore <- setdiff(tolower(ignore), bases)
  
  result <- sapply(
    apex::getSequences(as.multidna(x), simplify = FALSE),
    function(dna) {
      dna.mat <- tolower(as.character(as.matrix(dna)))
      site.freqs <- apply(dna.mat, 2, function(site) {
        site <- site[!site %in% ignore]
        table(factor(site, levels = bases))
      })
      colnames(site.freqs) <- 1:ncol(site.freqs)
      base.freqs <- table(factor(as.vector(dna.mat), levels = bases))
      ind.freqs <- t(sapply(as.character(as.list(dna)), function(x) {
        table(factor(tolower(x), levels = bases))
      }))
      list(
        site.freqs = site.freqs, 
        base.freqs = base.freqs, 
        ind.freqs = ind.freqs
      )
    },
    simplify = FALSE
  )
  
  if(length(result) == 1 & simplify) result[[1]] else result
}