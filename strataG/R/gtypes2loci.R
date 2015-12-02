#' @name gtypes2loci
#' @title Convert Between \code{gtypes} And \code{loci} objects.
#' @description Convert a \code{gtypes} object to a \code{\link[pegas]{loci}} object.
#' 
#' @param x a \linkS4class{gtypes} object.
#'  
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{initialize.gtypes}, \link{df2gtypes}, \link{sequence2gtypes}, 
#'   \link{gtypes2df}, \link{gtypes2genind}
#' 
#' @importFrom pegas as.loci
#' @export
#' 
gtypes2loci <- function(x) {
  mat <- as.matrix(x, one.col = TRUE, sep = "/")
  mat <- data.frame(cbind(as.character(strata(x)), mat))
  mat <- data.frame(lapply(mat, factor))
  as.loci(mat, col.pop = 1)
}
  