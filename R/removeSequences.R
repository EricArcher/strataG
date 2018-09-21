#' @title Remove Sequences
#' @description Remove sequences not used by samples.
#'   
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return a new \linkS4class{gtypes} object with unused sequences removed.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom methods new
#' @importFrom stats na.omit
#' @export
#' 
removeSequences <- function(g) {
  if(is.null(sequences(g))) return(g)
  haps <- getAlleleNames(g)
  g@sequences <- g %>% 
    getLocusNames() %>% 
    purrr::map(function(x) {
      sequences(g)[[x]][haps[[x]]]
    }) %>% 
    setNames(names(haps)) %>% 
    as.multidna()
  g
}
  
  