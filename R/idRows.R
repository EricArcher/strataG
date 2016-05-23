#' @title \code{gtypes} ID Rows
#' @description Return lookup matrix of row numbers or row names of 
#'   \code{@@loci} slot for sample ids.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param ids character vector of sample ids. If \code{NULL}, returns marix of 
#'   all ids.
#' @param as.names \code{FALSE} returns numeric matrix of row numbers. \code{TRUE} 
#'   returns character matrix of rownames.
#' 
#' @return a matrix of rownames or row numbers for sample ids in the \code{loci} 
#'   slot of a gtypes object.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#'
#' # lookup table of all samples by number
#' by.num <- idRows(msats.g)
#' head(by.num)
#' 
#' # lookup table of all samples by name
#' by.name <- idRows(msats.g, as.name = TRUE)
#' head(by.name)
#' 
#' # lookup table of selected ids
#' ran.ids <- sample(indNames(msats.g), 5)
#' ran.by.num <- idRows(msats.g, ids = ran.ids, as.name = TRUE)
#' ran.ids
#' 
#' @export
#'
idRows <- function(g, ids = NULL, as.names = FALSE) {
  g.ids <- indNames(g)
  if(!is.null(ids)) {
    missing <- setdiff(ids, g.ids)
    if(length(missing) > 0) {
      missing <- paste(missing, collapse = ", ")
      stop(paste("the following ids could not be found:", missing))
    }
  }
  id.vec <- sub("\\.[[:digit:]]*$", "", rownames(g@loci))
  mat.vec <- if(as.names) rownames(g@loci) else 1:nrow(g@loci)
  mat <- do.call(rbind, tapply(mat.vec, id.vec, function(i) i, simplify = FALSE))
  ids <- if(!is.null(ids)) unique(ids) else g.ids
  mat[ids, , drop = FALSE]
}