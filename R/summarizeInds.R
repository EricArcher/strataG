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
  g |> 
    zygosity() |> 
    dplyr::group_by(.data$id) |> 
    dplyr::summarize(
      num.loci.genotyped = sum(!is.na(.data$zyg)),
      num.loci.missing.genotypes = sum(is.na(.data$zyg)),
      pct.loci.missing.genotypes = mean(is.na(.data$zyg)),
      pct.loci.homozygous = sum(.data$zyg == 'hom', na.rm = TRUE) / 
        .data$num.loci.genotyped,
      .groups = 'drop'
    ) |> 
    dplyr::mutate(stratum = getStrata(g)[.data$id]) |> 
    dplyr::select('id', 'stratum', dplyr::everything()) |> 
    dplyr::arrange(.data$id, .data$stratum) |> 
    as.data.frame()
}
