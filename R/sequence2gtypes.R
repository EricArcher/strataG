#' @title Convert Sequences To \code{gtypes}
#' @description Create a \linkS4class{gtypes} object from sequence data.
#' 
#' @param x DNA sequences as a character matrix, a \code{DNAbin} object, 
#'   or \code{multidna} object.
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
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
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
  
  # get individual names
  ind.names <- unique(unlist(apex::getSequenceNames(x)))
  
  # check strata if provided
  if(!is.null(strata)) {
    # is vector named?
    if(!is.null(names(strata))) {
      # are all ids present?
      if(!all(ind.names %in% names(strata))) {
        extra <- setdiff(ind.names, names(strata))
        extra <- stats::setNames(rep(NA, length(extra)), extra)
        strata <- c(strata, extra)
      }
    } else {
      if(length(strata) != length(ind.names)) {
        warning("the length of 'strata' is not the same as the number of individuals. strata will be recycled.")
        strata <- rep(strata, length.out = length(ind.names))
      }
    }
  } else strata <- rep('Default', length(ind.names))
  
  # create and return gtypes object
  do.call(
    data.frame, 
    lapply(apex::getSequences(x, simplify = FALSE), function(dna) {
      x.labels <- labels(dna)
      x.labels[match(ind.names, x.labels)]
    })
  ) |> 
    stats::setNames(apex::getLocusNames(x)) |> 
    dplyr::mutate(
      id = ind.names, 
      strata = strata
    ) |> 
    dplyr::select(
      dplyr::all_of(c('id', 'strata')), 
      dplyr::everything()
    ) |> 
    df2gtypes( 
      ploidy = 1, 
      schemes = schemes, 
      sequences = x, 
      description = description, 
      other = if(is.null(other)) list() else other
    )
}
