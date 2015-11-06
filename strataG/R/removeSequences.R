#' @title Remove Sequences
#' @description Remove sequences not used by samples listed in \code{@@loci} 
#'   slot.
#'   
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return a new \linkS4class{gtypes} object with unused sequences removed.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom methods new
#' @export
#' 
removeSequences <- function(g) {
  if(is.null(sequences(g))) return(g)
  new.list <- lapply(locNames(g), function(x) {
    haps <- unique(as.character(loci(g)[[x]]))
    sequences(g, x)[haps, ]
  })
  names(new.list) <- locNames(g)
  g@sequences <- new("multidna", new.list)
  g
}
  
  