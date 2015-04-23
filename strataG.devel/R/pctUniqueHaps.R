#' @title Percent Unique Haplotypes
#' @description Calculate the percent of haplotypes that are unique.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
pctUniqueHaps <- function(g) { 
  if(g@ploidy != 1) stop("'g' is not haploid")
  apply(g@loci, 2, function(locus) sum(table(locus) == 1) / sum(!is.na(locus)))
}