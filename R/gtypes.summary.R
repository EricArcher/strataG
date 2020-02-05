#' @title Summarize gtypes Object
#' @description Generate a summary of a \code{gtypes} object.
#'  
#' @param object a \linkS4class{gtypes} object.
#' @param ... other arguments (ignored).
#' 
#' @return a list with the following elements:
#' \describe{
#'   \item{\code{num.ind}}{number of individuals}
#'   \item{\code{num.loc}}{number of loci}
#'   \item{\code{num.strata}}{number of strata}
#'   \item{\code{unstratified}}{number of unstratified samples}
#'   \item{\code{schemes}}{names of stratification schemes}
#'   \item{\code{allele.freqs}}{a list with tables of allele frequencies by strata}
#'   \item{\code{strata.smry}}{a by-strata data.frame summarizing haplotypes or loci}
#'   \item{\code{locus.smry}}{a data.frame summarizing each locus for 
#'     non-haploid objects, \code{NULL} for haploid objects}
#'   \item{\code{seq.smry}}{a summary of the sequence length and base frequencies}
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @name summary,gtypes-method
#' @aliases summary summary.gtypes
#' 
#' @export
#' 
setMethod("summary", "gtypes", function(object, ...) { 
  .printSmryHeader(list(
    description = getDescription(object),
    num.ind = getNumInd(object), 
    num.loc = getNumLoci(object), 
    num.strata = getNumStrata(object)
  )) 
  cat("\n")
  
  .summaryStats <- function(x) {
    apply(x, 2, function(z) {
      c(
        Min = min(z, na.rm = TRUE),
        Median = stats::median(z, na.rm = TRUE),
        Mean = mean(z, na.rm = TRUE),
        Max = max(z, na.rm = TRUE),
        NAs = sum(is.na(z))
      )
    })
  }
  
  print(cbind(
    summarizeInds(object) %>% 
      dplyr::select(-.data$id, -.data$stratum) %>% 
      .summaryStats(),
    summarizeLoci(object) %>% 
      dplyr::select(-.data$locus) %>% 
      .summaryStats()
  ))
})