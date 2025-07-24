#' @title Identify zygosity across loci
#' @description Identifies genotypes as heterozygous or homozygous.
#'   
#' @param x a \linkS4class{gtypes} object.
#' 
#' @return data frame that identifies if the genotype for an individual
#'   (\code{id}) at a locus (\code{locus}) is a heterozygote (\code{het}), 
#'   homozygote (\code{hom}) or missing (\code{NA}). Note that if haploid, 
#'   all non-missing values will be \code{hom}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' 
#' msats.zyg <- zygosity(msats.g)
#' head(msats.zyg, 10)
#' 
#' @export
#' 
zygosity <- function(x) {
  x@data |>
    dplyr::group_by(id, locus) |> 
    dplyr::summarize(
      zyg = ifelse(
        any(is.na(allele)),
        NA,
        ifelse(
          dplyr::n_distinct(allele, na.rm = TRUE) > 1,
          'het',
          'hom'
        )
      ),
      .groups = 'drop'
    )
}
