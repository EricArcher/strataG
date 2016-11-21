#' @title Stratify gtypes
#' @description Choose a new stratification scheme from the \code{schemes}
#'   slot in a \linkS4class{gtypes} object.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param scheme either the column name of a stratification scheme stored 
#'   in the data.frame of the \code{schemes} slot of \code{g}, or a vector or 
#'   factor identifying which stratum each sample belongs to.
#' @param drop remove samples not assigned to a stratum? (those assigned \code{NA} 
#'   in stratification scheme)
#'
#' @note If \code{scheme} is a vector or factor and has names, then the 
#'   they will be used to match with \code{\link{indNames}} of \code{g}. 
#'   Otherwise \code{scheme} should be the same length as the number of 
#'   samples in \code{g} or values in \code{scheme} will be recycled as 
#'   necessary.
#'   
#' @return A new \linkS4class{gtypes} object with an updated \code{strata} slot.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{schemes}
#'
#' @examples
#' data(msats.g)
#' msats.g
#' 
#' broad.msats <- stratify(msats.g, "broad")
#' broad.msats
#' 
#' @export
#' 
stratify <- function(g, scheme = NULL, drop = TRUE) {
  ids <- indNames(g)
  
  scheme <- if(is.null(scheme)) {
    rep("Default", length(ids))
  } else if(!(is.vector(scheme) | is.factor(scheme))) {
    stop("'scheme' must be a vector or a factor")
  } else if(length(scheme) == 1) {
    if(!scheme %in% colnames(schemes(g))) {
      stop(paste("scheme '", scheme, "' cannot be found", sep = ""))
    }
    g@schemes[ids, scheme]
  } else {
    if(length(scheme) != length(ids)) {
      warning(
        "'scheme' is not the same length as the number of samples. ",
        "values will be recycled."
      )
    }
    if(!is.null(names(scheme))) {
      scheme[ids]
    } else {
      rep(scheme, length.out = length(ids))
    }
  }
  scheme <- data.table(ids = ids, strata = as.character(scheme))
  g@data <- merge(scheme, g@data[, .SD, .SDcols = !"strata"], by = "ids", all.y = TRUE)
  
  if(drop) {
    g@data <- g@data[!is.na(strata)]
    g <- removeSequences(g)
  }
  g
}