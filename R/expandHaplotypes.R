#' @title Expand haplotypes
#' @description Expand haplotypes to a single sequence per individual.
#'
#' @param g a haploid \linkS4class{gtypes} object with sequences.
#'
#' @return a \code{gtypes} object with sequences expanded and renamed so there is
#'   one sequence per individual. Sequence names are set to individual sample IDs.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(dloop.g)
#' 
#' # Haplotypes have already been labelled
#' dloop.g
#' 
#' # Haplotypes expanded to individual sequences (num.alleles == num.samples)
#' expanded.g <- expandHaplotypes(dloop.g)
#' expanded.g
#' 
#' @export
#' 
expandHaplotypes <- function(g) {
  if(ploidy(g) != 1) stop("'g' must be a haploid gtypes object")
  dna.seqs <- sequences(g, as.haplotypes = FALSE)
  if(is.null(dna.seqs)) stop("'g' must have associated sequences")
  gen.data <- as.data.frame(g)
  for(x in getLocusNames(dna.seqs)) {
    to.replace <- which(!is.na(gen.data[, x]))
    gen.data[to.replace, x] <- as.character(gen.data[to.replace, 1])
  }
  
  new.g <- df2gtypes(
    gen.data, ploidy = 1, schemes = schemes(g), 
    sequences = dna.seqs, description = description(g),
    other = other(g)
  )
}



