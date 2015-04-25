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
removeSequences <- function(g) {
  new.list <- lapply(colnames(g@loci), function(x) {
    haps <- unique(as.character(g@loci[[x]]))
    g@sequences@dna[[x]][haps, ]
  })
  names(new.list) <- colnames(g@loci)
  g@sequences <- new("multidna", new.list)
  g
}
  
  