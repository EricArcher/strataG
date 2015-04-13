#' @title Stratify gtypes
#' @description Choose a new stratification scheme from the \code{schemes}
#'   slot in a \linkS4class{gtypes} object.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param scheme column name of a stratification scheme from the data.frame in
#'   the \code{schemes} slot of \code{g}.
#'
#' @return A new \linkS4class{gtypes} object with an updated \code{strata}
#'   slot.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export

stratify <- function(g, scheme = NULL) {
  ids <- indNames(g)
  scheme <- strata <- if(is.null(scheme)) {
    character(0)
  } else if(!(is.vector(scheme) | is.factor(scheme))) {
    stop("'scheme' must be a vector or a factor")
  } else if(length(scheme) == 1) {
    if(!scheme %in% colnames(g@schemes)) {
      stop(paste("scheme '", scheme, "' cannot be found", sep = ""))
    }
    g@schemes[ids, scheme]
  } else {
    if(length(scheme) != length(ids)) {
      warning(paste("'scheme' is not the same length as the number of samples.",
                    "values will be recycled"))
    }
    if(!is.null(names(scheme))) {
      factor(scheme[ids])
    } else {
      factor(scheme)
    }
  }
  names(scheme) <- ids
  g@strata <- scheme
  return(g)
}