#' @title Convert \code{gtypes} To \code{data.frame}
#' @description Create a data.frame from a \linkS4class{gtypes} object.
#'   
#' @param x a \linkS4class{gtypes} object.
#' @param one.col logical. If \code{TRUE}, then result has one column per 
#'   locus.
#' @param sep character to use to separate alleles if \code{one.col} is 
#'   \code{TRUE}.
#' @param ids logical. include a column for individual identifiers (\code{ids})?
#' @param strata logical. include a column for current statification (\code{strata})?
#' @param ... additional arguments ot be passed to or from methods.
#'   
#' @return \code{data.frame} with one row per sample.
#' 
#' @author Eric Archer \email{eric.archer@@noa.gov}
#' 
#' @seealso \link{df2gtypes}, \link{as.matrix.gtypes}
#' 
#' @aliases as.data.frame.gtypes
#' 
#' @export
#' 
setMethod("as.data.frame", "gtypes", function(x, one.col = FALSE, sep = "/", 
                                          ids = TRUE, strata = TRUE, ...) {
  as.data.frame(as.matrix(
    x, one.col = one.col, sep = sep, ids = ids, strata = strata
  ), ...)
})