#' @title Observed and Expected Heterozygosity 
#' @description Calculate observed heterozygosity for diploid data.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @note For a measure of haplotypic diversity (haploid "heterozygosity"), 
#'   use \code{exptdHet}. If \code{g} is a haploid object with sequences, make sure to run 
#'   \code{\link{labelHaplotypes}} if sequences aren't already grouped by 
#'   haplotype.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' 
#' # Expected heterozygosity
#' exptdHet(msats.g)
#' 
#' # Observed heterozygosity
#' obsvdHet(msats.g)
#' 
#' @name heterozygosity
#' 
NULL


#' @rdname heterozygosity
#' @importFrom swfscMisc diversity
#' @export
#' 
exptdHet <- function(g) {
  .applyPerLocus(swfscMisc::diversity, g)
}


#' @rdname heterozygosity
#' @importFrom stats na.omit
#' @export
#' 
obsvdHet <- function(g) {
  isHom <- function(x) {
    if(any(is.na(x))) NA else length(unique(x)) == 1
  }
  is.homzgt <- g@data[, lapply(.SD, isHom), .SDcols = !c("ids", "strata"), by = "ids"]
  is.homzgt[, 1 - sapply(.SD, mean, na.rm = TRUE), .SDcols = !"ids"]
}