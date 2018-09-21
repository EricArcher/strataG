#' @title Individual Summaries
#' @description Compile standard by-individual summaries.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return A data.frame with rows for each sample and columns containing:
#' \describe{
#'   \item{\code{id}}{The individual id}
#'   \item{\code{stratum}}{The stratum of the individual}
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
summarizeInd <- function(g) {
  g@data %>%
    dplyr::group_by(id, locus) %>% 
    dplyr::summarize(is.hmzgs = ifelse(
      any(is.na(allele)),
      NA,
      length(unique(allele)) == 1
    )) %>% 
    dplyr::group_by(id) %>% 
    dplyr::summarize(
      num.loci.genotyped = sum(!is.na(is.hmzgs)),
      num.loci.missing.genotypes = sum(is.na(is.hmzgs)),
      pct.loci.missing.genotypes = mean(is.na(is.hmzgs)),
      pct.loci.homozygous = sum(is.hmzgs, na.rm = TRUE) / dplyr::n()
    ) %>% 
    dplyr::mutate(stratum = strata(g)[id]) %>% 
    dplyr::select(id, stratum, dplyr::everything()) %>% 
    dplyr::arrange(stratum, id) %>% 
    as.data.frame()
}
