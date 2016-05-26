#' @title Convert Between \code{gtypes} And \code{phyDat} objects.
#' @description Convert a \code{gtypes} object to a \code{\link[phangorn]{phyDat}} object.
#' 
#' @param x a \linkS4class{gtypes} or \code{phyDat} formatted object.
#' @param locus name or number of locus to convert.
#' @param ... optional arguments passed to \code{\link{sequence2gtypes}}.
#'  
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{initialize.gtypes}, \link{df2gtypes}, 
#'   \link{sequence2gtypes}, \link{as.data.frame.gtypes}, 
#'   \link{gtypes2genind}, \link{gtypes2loci}
#' 
#' @examples
#' data(dloop.g)
#' 
#' # Convert to phDat
#' pd <- gtypes2phyDat(dloop.g)
#' pd
#' 
#' # Convert to gtypes
#' gt <- phyDat2gtypes(pd)
#' gt 
#' 
#' @name gtypes2phyDat 
#' @importFrom phangorn as.phyDat
#' @export
#' 
gtypes2phyDat <- function(x, locus = 1) {
  x <- getSequences(sequences(x), loci = locus, simplify = TRUE) 
  as.phyDat(as.matrix(x))
}

#' @rdname gtypes2phyDat
#' @export
#' 
phyDat2gtypes <- function(x, ...) {
  sequence2gtypes(as.DNAbin(x), ...)
}
  