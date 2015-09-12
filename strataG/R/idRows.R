#' @title \code{gtypes} ID Rows
#' @description Return row numbers of \code{@@loci} slot for a set of
#'   sample ids.
#'   
#' @param ids character vector of sample ids.
#' @param rowNames a vector of rownames from the \code{@@loci} slot of a 
#'   \linkS4class{gtypes} object.
#' 
#' @return a vector of row numbers in \code{g@@loci} where data for \code{ids}
#'   are stored.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.msats)
#' data(dolph.strata)
#' msats.merge <- merge(dolph.strata[, c("ids", "fine")], dolph.msats, all.y = TRUE)
#' msats <- df2gtypes(msats.merge, ploidy = 2)
#' 
#' ran.ids <- sample(indNames(msats), 5)
#' idRows(ran.ids, rownames(loci(msats)))
#' 
#' @export
#' 
idRows <- function(ids, rowNames) {
  ids <- as.character(ids)
  id.rows <- sub("\\.[[:digit:]]*$", "", rowNames)
  rows <- sapply(ids, function(x) which(id.rows == x))
  if(is.vector(rows)) {
    matrix(rows, ncol = 1, dimnames = list(ids, NULL))
  } else {
    t(rows)
  }
}