#' @title Convert to multidna
#' @description Convert a set of sequences to a multidna object if possible.
#' 
#' @param x a valid set of sequences: character matrix, list of 
#'   character vectors, \code{\link{DNAbin}} object or list of them,
#'   \linkS4class{gtypes} object, or \linkS4class{multidna} object.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom methods new
#' @export
#' 
as.multidna <- function(x) {
  # multidna
  if(inherits(x, "multidna")) return(x)
  
  # gtypes with sequences
  if(inherits(x, "gtypes")) {
    if(is.null(x@sequences)) {
      stop("the gtypes object does not contain sequences")
    }
    return(x@sequences)
  }
  
  # a DNAbin
  if(inherits(x, "DNAbin")) return(new("multidna", list(as.matrix(x))))
  
  # character matrix or list of character vectors
  if(is.character(x) | (is.list(x) & is.character(x[[1]]))) {
    return(new("multidna", list(as.matrix(as.DNAbin(x)))))
  }
  
  # list of DNAbin
  if(is.list(x) & inherits(x[[1]], "DNAbin")) {
    x <- lapply(x, as.matrix)
    return(new("multidna", x))
  }
  
  stop("'x' must be a valid set of sequences")
}