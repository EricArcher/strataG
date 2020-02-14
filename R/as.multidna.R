#' @title Convert to multidna
#' @description Convert a set of sequences to a multidna object if possible.
#' 
#' @param x a valid set of sequences: character matrix, list of 
#'   character vectors, \code{\link{DNAbin}} object or list of them,
#'   \linkS4class{gtypes} object, or \linkS4class{multidna} object.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{getSequences}
#' 
#' @examples 
#' # convert list of character vectors
#' data(dolph.seqs)
#' list.mdna <- as.multidna(dolph.seqs)
#' list.mdna
#' 
#' # convert gtypes object
#' data(dloop.g)
#' gtype.mdna <- as.multidna(dloop.g)
#' gtype.mdna
#' 
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
    return(getSequences(x, as.multidna = TRUE))
  }
  
  # a DNAbin
  if(inherits(x, "DNAbin")) return(methods::new("multidna", list(as.matrix(x))))
  
  # character matrix or list of character vectors
  if(is.character(x) | (is.list(x) & all(sapply(x, is.character)))) {
    x <- list(ape::as.DNAbin(x))
  }
  
  # list of DNAbin
  if(is.list(x) & all(sapply(x, function(elem) inherits(elem, "DNAbin")))) {
    x <- sapply(x, as.matrix, simplify = FALSE)
    if(is.null(names(x))) names(x) <- paste("gene", 1:length(x), sep = "")
    return(methods::new("multidna", x))
  }
  
  stop("'x' must be a valid set of sequences")
}