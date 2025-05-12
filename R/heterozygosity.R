#' @title Heterozygosity 
#' @description Calculate observed and heterozygosity.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results by strata?
#' @param type return \code{expected} or \code{observed} heterozygosity
#' 
#' @note If \code{g} is a haploid object with sequences, the value for 
#'   expected heterozygosity (= haplotpyic diversity) will be returned.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' 
#' # Expected heterozygosity
#' heterozygosity(msats.g, type = "expected")
#' 
#' # Observed heterozygosity by strata
#' heterozygosity(msats.g, FALSE, "observed")
#' 
#' @export
#' 
heterozygosity <- function(g, by.strata = FALSE, type = c("expected", "observed")) {
  if(getPloidy(g) == 1) type <- "expected"
  g <- .checkHapsLabelled(g)
  
  result <- switch(
    match.arg(type),
    expected = .applyPerLocus(sprex::diversity, g, by.strata = by.strata, type = "unb.gini") |> 
      dplyr::rename(exptd.het = .data$value),
    observed = {
      is.het <- if(by.strata) {
        g@data |> 
          dplyr::group_by(.data$stratum, .data$locus, .data$id) |> 
          dplyr::summarize(
            is.het = dplyr::n_distinct(.data$allele) > 1, 
            .groups = "drop_last"
          )
      } else {        
        g@data |> 
          dplyr::group_by(.data$locus, .data$id) |> 
          dplyr::summarize(
            is.het = dplyr::n_distinct(.data$allele) > 1, 
            .groups = "drop_last"
          ) 
      }
      dplyr::summarize(
        is.het,
        obsvd.het = mean(.data$is.het, na.rm = TRUE),
        .groups = "drop"
      )
    }
  )
  
  result <- if("stratum" %in% colnames(result)) {
    dplyr::arrange(result, .data$stratum, .data$locus)
  } else {
    dplyr::arrange(result, .data$locus)
  }
  
  as.data.frame(result)
}
