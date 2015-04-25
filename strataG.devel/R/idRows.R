#' @title ID Rows
#' @description Return row numbers of \code{@@loci} slot for a set of
#'   sample ids.
#'   
#' @param ids character vector of sample ids.
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return a vector of row numbers in \code{g@@loci} where data for \code{ids}
#'   are stored.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
idRows <- function(ids, g) {
  ids <- as.character(ids)
  id.rows <- sub("\\.[[:digit:]]*$", "", rownames(g@loci))
  unlist(sapply(ids, function(x) which(id.rows == x)))
}