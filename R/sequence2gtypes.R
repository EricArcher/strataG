#' @title Convert Sequences To \code{gtypes}
#' @description Create a \linkS4class{gtypes} object from sequence data.
#' 
#' @param x DNA sequences as a character matrix, a \code{\link{DNAbin}} object, 
#'   or \linkS4class{multidna} object.
#' @param strata a vector or factor giving stratification for each sequence. If 
#'   not provided all individuals are assigned to the same stratum (Default).
#' @param seq.names names for each set of sequences. If not provided default names 
#'   are generated.
#' @param schemes an optional data.frame of stratification schemes.
#' @param description an optional label for the object.
#' @param other a list to carry other related information (optional).
#' 
#' @return a \linkS4class{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{gtypes.initialize}, \link{as.matrix.gtypes}, 
#'   \link{as.data.frame.gtypes}, \link{gtypes2genind}, \link{gtypes2loci},
#'   \link{gtypes2phyDat}
#' 
#' @examples
#' #--- create a haploid sequence (mtDNA) gtypes object
#' data(dolph.strata)
#' data(dolph.seqs)
#' strata <- dolph.strata$fine
#' names(strata) <- dolph.strata$ids
#' dloop.fine <- sequence2gtypes(dolph.seqs, strata, seq.names = "dLoop", 
#' description = "dLoop: fine-scale stratification")
#' 
#' @export
#' 
sequence2gtypes <- function(x, strata = NULL, seq.names = NULL, schemes = NULL,
                            description = NULL, other = NULL) {
  # convert sequences
  x <- as.multidna(x)
  
  # check seq.names
  if(!is.null(seq.names)) {
    if(length(seq.names) != getNumLoci(x)) {
      stop("length of 'seq.names' is not equal to number of genes")
    }
    x <- apex::"setLocusNames<-"(x, seq.names)
  }
  
  # create gen.data data.frame
  ind.names <- unique(unlist(apex::getSequenceNames(x)))
  gen.data <- do.call(
    data.frame, 
    lapply(apex::getSequences(x, simplify = FALSE), function(dna) {
      x.labels <- labels(dna)
      x.labels[match(ind.names, x.labels)]
    })
  )
  colnames(gen.data) <- apex::getLocusNames(x)
  rownames(gen.data) <- ind.names
  
  # return new gtypes object
  methods::new(
    "gtypes", 
    gen.data = gen.data, 
    ploidy = 1, 
    strata = strata,
    schemes = schemes, 
    sequences = x, 
    description = description, 
    other = if(is.null(other)) list() else other
  )
}