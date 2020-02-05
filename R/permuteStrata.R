#' @title Permute strata
#' @description Permute the strata slot within a \linkS4class{gtypes} object.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return a \linkS4class{gtypes} object with the strata randomly permuted.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#' summary(msats.g)
#' 
#' ran.msats <- permuteStrata(msats.g)
#' summary(ran.msats)
#'
#' @export
#' 
permuteStrata <- function(g) {
  st <- getStrata(g)
  no.nas <- st[!is.na(st)]
  no.na.perm <- sample(no.nas)
  names(no.na.perm) <- names(no.nas)
  setStrata(g) <- stats::setNames(no.na.perm[names(st)], names(st))
  g
}