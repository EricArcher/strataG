#' @title Locus Summaries
#' @description Compile standard by-locus summaries.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical. If \code{TRUE}, return a list of summary matrices 
#'   for each stratum.
#' @param ... arguments to be passed on to summary functions.
#' 
#' @return A matrix with rows for each locus and columns containing summary
#'   statistics.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.msats)
#' data(dolph.strata)
#' msats.merge <- merge(dolph.strata[, c("ids", "fine")], dolph.msats, 
#'   all.y = TRUE)
#' msats.g <- df2gtypes(msats.merge, ploidy = 2)
#' summarizeLoci(msats.g)
#' 
#' @export
#' 
summarizeLoci <- function(g, by.strata = FALSE, ...) {
  summary.stats <- function(x) {
    cbind(
      num.missing = numMissing(x),
      num.alleles = numAlleles(x),
      allelic.richness = allelicRichness(x),
      expected.het = exptdHet(x),
      observed.het = obsvdHet(x)
    )
  }
  
  if(by.strata) {
    lapply(strataSplit(g), summary.stats)
  } else summary.stats(g)
}