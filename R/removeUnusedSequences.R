#' @title Remove Unused Sequences
#' @description Remove sequences not used by samples.
#'   
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return a new \linkS4class{gtypes} object with unused sequences removed.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
removeUnusedSequences <- function(g) {
  if(is.null(getSequences(g))) return(g)
  haps <- getAlleleNames(g)
  g@sequences <- g |> 
    getLociNames() |> 
    purrr::map(function(x) {
      getSequences(g, as.haplotypes = TRUE, seqName = x, simplify = TRUE)[haps[[x]]]
    }) |> 
    stats::setNames(names(haps)) |> 
    as.multidna()
  g
}
  
  
