#' @name as.array.gtypes
#' @title Convert \code{gtypes} To \code{array}
#' @description Create a 3-dimensional array from a \linkS4class{gtypes} object.
#'   
#' @param x a \linkS4class{gtypes} object.
#' @param ids vector of individual ids.
#' @param loci vector of loci.
#' @param drop if \code{TRUE} the array is coerced to the lowest possible dimension.
#' @param ... additional arguments to be passed to or from methods.
#'   
#' @return A three dimensional \code{array} with with dimeinsions of: 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link[strataG]{as.data.frame} \link[strataG]{as.matrix}
#' 
#' @examples 
#' data(msats.g)
#' msats.arr <- as.array(msats.g)
#' 
#' msats.arr
#' 
#' @aliases as.array,gtypes-method as.array.gtypes as.array
#' @importFrom methods setMethod
#'
#' @export
#' 
setMethod(
  "as.array", "gtypes", 
  function(x, drop = TRUE, ...) {
  
  # create 3-D array
  arr <- array(
    dim = c(length(ids), length(loci), ploidy(x)),
    dimnames = list(ids = ids, loci = loci, NULL)
  )
  for(i in ids) arr[i, , ] <- t(as.matrix(x@data[i, loci, with = FALSE]))
  if(drop) drop(arr) else arr
})