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
#' summarizeInds(msats.g)
#' 
#' @export
#' 
summarizeInds <- function(g) {
  g@data %>%
    dplyr::group_by(.data$id, .data$locus) %>% 
    dplyr::summarize(is.hmzgs = ifelse(
      any(is.na(.data$allele)),
      NA,
      length(unique(.data$allele)) == 1
    )) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(.data$id) %>% 
    dplyr::summarize(
      num.loci.genotyped = sum(!is.na(.data$is.hmzgs)),
      num.loci.missing.genotypes = sum(is.na(.data$is.hmzgs)),
      pct.loci.missing.genotypes = mean(is.na(.data$is.hmzgs)),
      pct.loci.homozygous = sum(.data$is.hmzgs, na.rm = TRUE) / dplyr::n()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(stratum = getStrata(g)[.data$id]) %>% 
    dplyr::select(.data$id, .data$stratum, dplyr::everything()) %>% 
    dplyr::arrange(.data$stratum, .data$id) %>% 
    as.data.frame()
}
