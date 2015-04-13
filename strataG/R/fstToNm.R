#' @title Nm from Fst
#' @description Calculate Nm for a given value of Fst.
#'
#' @param fst estimate of Fst between populations.
#' @param diploid logical. Is Fst from diploid (TRUE) or haploid (FALSE) data?
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export

fstToNm <- function(fst, diploid = TRUE) {
  ((1 / fst) - 1) / ifelse(diploid, 4, 2)
}