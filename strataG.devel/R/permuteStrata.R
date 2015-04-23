#' @title Permute strata
#' @description Permute the strata slot within a \linkS4class{gtypes} object.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.msats)
#' data(dolph.strata)
#' msats.merge <- merge(dolph.strata[, c("ids", "fine")], dolph.msats, 
#'   all.y = TRUE, description = date())
#' msats <- df2gtypes(msats.merge, ploidy = 2)
#' summary(msats)
#' ran.msats <- permuteStrata(msats)
#' summary(ran.msats)
#'
#' @export

permuteStrata <- function(g) {
  strata <- g@strata
  no.nas <- strata[!is.na(strata)]
  no.na.sample <- sample(as.character(no.nas))
  names(no.na.sample) <- names(no.nas)
  g@strata <- factor(no.na.sample[names(strata)])
  g
}