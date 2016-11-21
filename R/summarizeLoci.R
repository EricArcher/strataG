#' @title Locus Summaries
#' @description Compile standard by-locus summaries.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical. If \code{TRUE}, return a list of summary matrices 
#'   for each stratum.
#' @param ... arguments to be passed on to summary functions.
#' 
#' @return A matrix with rows for each locus and columns containing summaries of:
#' \describe{
#'   \item{\code{num.genotyped}}{The number of samples genotyped}
#'   \item{\code{prop.genotyped}}{The proportion of samples genotyped}
#'   \item{\code{num.alleles}}{The number of alleles in the locus}
#'   \item{\code{allelic.richness}}{The allelic richness of the locus}
#'   \item{\code{prop.unique.alleles}}{Proportion of alleles found in a single sample}
#'   \item{\code{expt.heterozygosity}}{Expected heterozygosity}
#'   \item{\code{obsvd.heterozygosity}}{Observed heterozygosity}
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#' 
#' summarizeLoci(msats.g)
#' 
#' @export
#' 
summarizeLoci <- function(g, by.strata = FALSE, ...) {
  smry.func <- function(x) {
    n.gtyped <- numGenotyped(x)
    cbind(
      num.genotyped = n.gtyped,
      prop.genotyped = n.gtyped / nInd(x),
      num.alleles = numAlleles(x),
      allelic.richness = allelicRichness(x),
      prop.unique.alleles = propUniqueAlleles(x),
      exptd.heterozygosity = exptdHet(x),
      obsvd.heterozygosity = obsvdHet(x)
    )
  }
  
  if(by.strata) {
    lapply(strataSplit(g), smry.func)
  } else smry.func(g)
}