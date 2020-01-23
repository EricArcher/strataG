#' @title Nucleotide Diversity
#' @description Calculate nucleotide diversity for set of haplotypes.
#' 
#' @param x a set of sequences or a \linkS4class{gtypes} object with sequences.
#' @param bases nucleotides to consider when calculating diversity.
#' @param simplify if \code{TRUE} and only one loci exists, return a vector, 
#'   otherwise, a list of vectors with one element per locus will be returned.
#' 
#' @return Vector of nucleotide diversity by site. 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dloop.g)
#' 
#' nd <- nucleotideDiversity(dloop.g)
#' quantile(nd)
#' 
#' @export
#' 
nucleotideDiversity <- function(x, bases = c("a", "c", "g", "t"), simplify = TRUE) {
  bases <- tolower(bases)
  
  x <- if(is.gtypes(x)) {
    getSequences(x, as.haplotypes = FALSE, as.multidna = TRUE)
  } else {
    as.multidna(x)
  }
  
  result <- sapply(
    apex::getSequences(x, simplify = FALSE), 
    function(dna) {  
      site.div <- dna %>% 
        as.matrix() %>% 
        as.character() %>% 
        apply(2, function(b) swfscMisc::diversity(b[b %in% bases]))
      stats::setNames(site.div, 1:length(site.div))
    },
    simplify = FALSE
  )
  
  if(length(result) == 1 & simplify) result[[1]] else result
}