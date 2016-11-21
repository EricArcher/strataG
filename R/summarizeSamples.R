#' @title Sample Summaries
#' @description Compile standard by-sample summaries.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param sort.by.strata logical. Sort data.frame by strata?
#' 
#' @return A data.frame with rows for each sample and columns containing:
#' \describe{
#'   \item{\code{id}}{The sample id}
#'   \item{\code{strata}}{The stratum of the sample}
#'   \item{\code{num.loci.missing.genotypes}}{The number of genotypes missing}
#'   \item{\code{pct.loci.missing.genotypes}}{The proportion of genotypes missing}
#'   \item{\code{pct.loci.homozygous}}{The proportion of loci homozygous}
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' 
#' summarizeSamples(msats.g)
#' 
#' @export
#' 
summarizeSamples <- function(g, sort.by.strata = FALSE) {
  result <- do.call(rbind, lapply(indNames(g), function(id) {
    smry <- sapply(locNames(g), function(loc) {
      gt <- as.array(g, id, loc)
      missing <- any(is.na(gt))
      hmzgs <- if(missing) NA else length(unique(unlist(gt))) == 1
      c(missing = missing, hmzgs = hmzgs)
    })
    
    missing <- sum(smry["missing", ], na.rm = TRUE)
    result <- data.frame(
      id = id, 
      strata = strata(g)[id],
      num.loci.missing.genotypes = missing,
      pct.loci.missing.genotypes = missing / nLoc(g)
    )
    result$pct.loci.homozygous <- if(ploidy(g) > 1) {
      mean(smry["hmzgs", ], na.rm = TRUE)
    } else NULL
    result
  }))
  rownames(result) <- indNames(g)
  result$id <- as.character(result$id)
  
  if(sort.by.strata) result <- result[order(result$strata, result$id), ]
  result
}
  