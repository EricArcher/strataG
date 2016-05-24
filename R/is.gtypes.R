#' @title Test if object is \code{gtypes}
#' @description Test if object is \code{gtypes}
#' 
#' @param x R object to be tested.
#' 
#' @return Logical stating if 'x' is a \code{\linkS4class{gtypes}} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' is.gtypes(msats.g) # TRUE
#' 
#' data(dolph.msats)
#' is.gtypes(dolph.msats) # FALSE
#' 
#' @export
#' 
is.gtypes <- function(x) inherits(x, "gtypes")