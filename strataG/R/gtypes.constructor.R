#' @aliases initialize,gtypes-methods new.gtypes
#'
#' @param .Object
#' @param loci
#' @param ploidy
#' @param ind.names
#' @param schemes
#' @param sequences
#' @param seq.names
#' @param description
#' @param other
#' @param strata
#'
#' @seealso the \linkS4class{gtypes} class
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'


setMethod("initialize", "gtypes",
          function(.Object, loci, ploidy, ind.names = NULL,
                   schemes = NULL, sequences = NULL, seq.names = NULL,
                   description = NULL, other = NULL, strata = NULL) {
  # check loci
  if(is.vector(loci)) loci <- cbind(loci)
  if(!(is.matrix(loci) | is.data.frame(loci))) {
    stop("'loci' is not a vector, matrix, or data.frame")
  }

  # check ploidy
  if(!identical(trunc(ploidy), ploidy)) {
    stop("'ploidy' must be an integer value.")
  }
  if(ncol(loci) %% ploidy != 0) {
    stop("the number of columns in 'loci' is not a multiple of 'ploidy'.")
  }
  ploidy <- as.integer(ploidy)

  # check ind.names
  ind.names <- if(!is.null(ind.names)) {
    if(!is.vector(ind.names)) stop("'ind.names' must be a vector")
    if(length(ind.names) != nrow(loci)) {
      stop("the length of 'ind.names' must equal the number of rows in 'loci'")
    }
    as.character(ind.names)
  } else {
    if(!is.null(rownames(loci))) rownames(loci) else 1:nrow(loci)
  }
  rownames(loci) <- ind.names

  # check schemes
  if(is.null(schemes)) {
    schemes <- data.frame()
  } else {
    if(!(is.matrix(schemes) | is.data.frame(schemes))) {
      stop("'schemes' is not a matrix or data.frame")
    } else {
      schemes <- as.data.frame(schemes)
      for(i in 1:ncol(schemes)) schemes[, i] <- factor(schemes[, i])
      if(is.null(rownames(schemes))) {
        if(nrow(schemes) == length(ind.names)) {
          rownames(schemes) <- ind.names
        } else {
          stop("'schemes' doesn't have rownames and not as long as 'ind.names'")
        }
      } else {
        if(length(intersect(rownames(schemes), rownames(loci))) == 0) {
          stop("no rownames in 'schemes' are in 'loci'")
        }
      }
    }
  }

  # check description
  if(is.null(description)) {
    description <- paste("gtypes object created on", date())
  }

  # format loci
  locus.cols <- matrix(1:ncol(loci), nrow = ploidy)
  new.loci <- do.call(data.frame, lapply(1:ncol(locus.cols), function(i) {
    factor(unlist(loci[, locus.cols[, i]], use.names = FALSE))
  }))
  colnames(new.loci) <- if(is.null(colnames(loci))) {
    # return generic names if no colnames assigned
    locus.nums <- formatC(1:ncol(new.loci),
                          digits = floor(log10(n.loc)),
                          flag = "0"
    )
    paste("Locus", locus.nums, sep = "_")
  } else parseLocusNames(colnames(loci), ploidy)
  rownames(new.loci) <- paste(rep(rownames(loci), ploidy),
                              rep(1:ploidy, each = nrow(loci)), sep = ".")
  loci <- new.loci

  # check sequences
  sequences <- if(!is.null(sequences)) {
    switch(class(sequences)[1],
           DNAbin = new("multidna", list(sequences)),
           multidna = sequences,
           new("multidna", list(as.DNAbin(sequences)))
    )
  } else {
    new("multidna")
  }

  # compare sequences to loci
  if(length(sequences@dna) > 0) {
    if(length(sequences) != ncol(loci)) {
      stop("the number of sets of sequences is not equal to the number of loci")
    }
    if(is.null(seq.names)) {
      seq.names <- paste("Locus", 1:ncol(loci), sep = "_")
    } else if(length(seq.names) != ncol(loci)) {
      stop("the length of 'seq.names' is not equal to the number of loci")
    }
    colnames(loci) <- names(sequences@dna) <- seq.names
  }

  # create gtypes object
  g <- new("gtypes", loci = loci, ploidy = ploidy, sequences = sequences,
           schemes = schemes, strata = factor(),
           description = description, other = other)

  if(!is.null(strata)) g <- stratify(g, strata)
  return(g)
}
