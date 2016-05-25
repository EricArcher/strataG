#' @title Observed and Expected Heterozygosity 
#' @description Calculate observed heterozygosity for diploid data.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @note For a measure of haplotypic diversity (haploid "heterozygosity"), 
#'   use \code{exptdHet}.
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
  apply(loci(g), 2, swfscMisc::diversity)
}


#' @rdname heterozygosity
#' @importFrom stats na.omit
#' @export
#' 
obsvdHet <- function(g) {
  apply(loci(g), 2, function(locus) {
    locus <- na.omit(matrix(locus, ncol = ploidy(g)))
    is.homozgt <- apply(locus, 1, function(x) length(unique(x)) == 1)
    1 - (sum(is.homozgt) / nrow(locus))
  })
}