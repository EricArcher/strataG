#' @title Show a gtypes object
#' @description Show a gtypes object
#' 
#' @param object a \linkS4class{gtypes} object.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @importFrom methods setMethod
#' 
setMethod("show", "gtypes", function(object) {
  x <- object
  print(summary(x))
})
  