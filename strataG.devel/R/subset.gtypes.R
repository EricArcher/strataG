#' @title Subset gtypes
#' @description Return subsets of a gtypes object based on specified
#'   strata, individuals, or loci.
#'   
#' @param x a \linkS4class{gtypes} object.
#' @param strata a character vector designating individuals of which strata 
#'   to select.
#' @param ids a character vector designating which individuals to select.
#' @param loci a character vector designating which loci to select.
#' @param remove.sequences logical. If \code{TRUE} any sequences not referenced 
#'   in selected samples will not be in the returned object.
#' @param ... optional arguments passed through generic (ignored).
#' 
#' @return a \linkS4class{gtypes} object containing only individuals and loci
#'   that match the above criteria. If no individuals match (e.g., all ids are
#'   in non-selected strata), the return is \code{NULL}.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @aliases subset,gtypes
#' @export
#' 
setMethod("subset", "gtypes",  
  function(x, strata = NULL, ids = NULL, loci = NULL, 
           remove.sequences = FALSE, ...) {
    # check that ids can be found
    x.ids <- indNames(x)
    if(!is.null(ids)) {
      if(!(is.atomic(ids) & is.character(ids))) {
        stop("'ids' must be a character vector.")
      }
      ids <- unique(ids)
      ids.found <- ids %in% x.ids
      if(!all(ids.found)) {
        missing.ids <- paste(ids[!(ids.found)], collapse = ", ")
        stop("The following ids could not be found:", missing.ids)
      }
    } else ids <- x.ids
    
    # check that strata can be found
    if(!is.null(strata)) {
      if(!(is.atomic(strata) & is.character(strata))) {
        stop("'strata' must be a character vector.")
      }
      strata <- unique(strata)
      strata.found <- strata %in% strataNames(x)
      if(!all(strata.found)) {
        missing.strata <- paste(strata[!(strata.found)], collapse = ", ")
        stop("The following strata could not be found:", missing.strata)
      }
      ids <- intersect(ids, x.ids[strata(x) %in% strata])
    }
    
    # check that loci can be found
    x.loci <- locNames(x)
    loci <- if(!is.null(loci)) {
      if(!(is.atomic(loci) & is.character(loci))) {
        stop("'loci' must be a character vector.")
      }
      loci.found <- loci %in% x.loci
      if(!all(loci.found)) {
        missing.loci <- paste(loci[!(loci.found)], collapse = ", ")
        stop(paste("The following loci could not be found:", missing.loci))
      }
      loci[loci %in% x.loci]
    } else x.loci
    
    x@loci <- x@loci[idRows(ids, x), loci, drop = FALSE]
    x@strata <- droplevels(x@strata[ids])
    if(remove.sequences) x <- removeSequences(x)
    x
})

