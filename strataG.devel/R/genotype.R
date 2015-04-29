#' @title Genotype
#' @description Get genotype of a specified set of individuals and loci
#' 
#' @param ids character vector of individual ids.
#' @param loci character vector of loci.
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return a data.frame of the genotypes requested.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
genotype <- function(ids, loci, g) {
  if(!all(ids %in% indNames(g))) stop("some 'ids' not found in g")
  if(!all(loci %in% locNames(g))) stop("some 'loci' not found in g")
  g@loci[idRows(ids, rownames(g@loci)), loci]
}
  