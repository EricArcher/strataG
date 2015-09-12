#' @title Convert Sequences To \code{gtypes}
#' @description Load sequence data into a \linkS4class{gtypes} object. 
#' 
#' @param x DNA sequences as a character matrix, a \code{\link{DNAbin}} object, 
#'   or \linkS4class{multidna} object.
#' @param strata a vector or factor giving stratification for each sequence.
#' @param seq.names names for each set of sequences.
#' @param description a label for the object (optional).
#' @param other a slot to carry other related information - unused in package
#'   analyses (optional).
#' 
#' @return a \linkS4class{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' #--- create a haploid sequence (mtDNA) gtypes object
#' data(dolph.strata)
#' data(dolph.seqs)
#' strata <- dolph.strata$fine
#' names(strata) <- dolph.strata$ids
#' dloop.fine <- sequence2gtypes(dolph.seqs, strata, seq.names = "dLoop",
#'   description = "dLoop: fine-scale stratification")
#' 
#' @importFrom methods new
#' @export
#' 
sequence2gtypes <- function(x, strata = NULL, seq.names = NULL, 
                            description = NULL, other = NULL) {
  # convert sequences
  x <- as.multidna(x)
  
  # check seq.names
  if(!is.null(seq.names)) {
    if(length(seq.names) != length(x@dna)) {
      stop("length of 'seq.names' is not equal to number of genes")
    }
    names(x@dna) <- seq.names
  }
  
  # create gen.data data.frame
  ind.names <- x@labels
  gen.data <- do.call(data.frame, lapply(x@dna, function(dna) {
    x.labels <- labels(dna)
    factor(x.labels[match(ind.names, x.labels)])
  }))
  colnames(gen.data) <- labels(x@dna)
  rownames(gen.data) <- ind.names
  
  # return new gtypes object
  new("gtypes", gen.data = gen.data, ploidy = 1, strata = strata,
      sequences = x, description = description, other = other
  )
}