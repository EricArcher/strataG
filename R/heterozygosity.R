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
    expected = .applyPerLocus(
      sprex::diversity, 
      g, 
      by.strata = by.strata, 
      type = "unb.gini"
    ) |> 
      dplyr::rename(exptd.het = .data$value),
    
    observed = {
      zyg.df <- zygosity(g)
      
      zyg.df <- if(by.strata) {
        zyg.df |> 
          dplyr::mutate(stratum = getStrata(g)[.data$id]) |> 
          dplyr::group_by(.data$stratum, .data$locus)
      } else dplyr::group_by(zyg.df, .data$locus)
      
      zyg.df |> 
        dplyr::summarize(
          obsvd.het = sum(!is.na(.data$zyg) & .data$zyg == 'het') / 
            sum(!is.na(.data$zyg)),
          .groups = 'drop'
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
