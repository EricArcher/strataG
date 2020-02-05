#' @title Allele Frequencies
#' @description Calculate allele frequencies or proportions for each locus.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical determining if results should be returned by strata?
#' @param type return counts (\code{"freq"}) or proportions (\code{"prop"})
#'
#' @return A list of allele frequencies for each locus. Each element is a
#'   vector (\code{by.strata = FALSE}) or matrix (\code{by.strata = TRUE}) 
#'   with the frequency or proportion of each allele. 
#'   
#' @note If \code{g} is a haploid object with sequences, the function will run 
#'   \code{\link{labelHaplotypes}} if sequences aren't already grouped by 
#'   haplotype. The \code{gtypes} object used with haplotype assignments and 
#'   unassigned individuals will be stored in \code{attr(*, "gtypes")}.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(msats.g)
#' 
#' f <- alleleFreqs(msats.g)
#' f$D11t # Frequencies for Locus D11t
#' 
#' f.pop <- alleleFreqs(msats.g, TRUE, "prop")
#' f.pop$EV94[, "Coastal"] # Proportions of EV94 alleles in the Coastal population
#' 
#' @export
#' 
alleleFreqs <- function(g, by.strata = FALSE, type = c("freq", "prop")) {
  g <- .checkHapsLabelled(g)
  
  # create list of allele frequencies for each locus
  af <- if(by.strata) {
    lapply(split(g@data, g@data$locus), function(x) table(x$allele, x$stratum))
  } else {
    lapply(split(g@data, g@data$locus), function(x) table(x$allele))
  }
  
  # if proportions requested divide frequencies by sum of frequencies (in each stratum)
  if(match.arg(type) == "prop") {
    af <- lapply(af, function(x) {
      if(length(dim(x)) == 1) x / sum(x) else t(t(x) / colSums(x))
    })
  }
  
  unassigned <- getOther(g, "haps.unassigned")
  if(!is.null(unassigned)) attr(af, "gtypes") <- g
  
  af
}