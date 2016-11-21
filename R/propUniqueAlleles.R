#' @title Proportion Unique Alleles
#' @description Calculate the proportion of alleles that are unique.
#' 
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return a vector of the proportion of unique (occuring only in one individual) 
#'   alleles for each locus.
#' 
#' @note If \code{g} is a haploid object with sequences, make sure to run 
#'   \code{\link{labelHaplotypes}} if sequences aren't already grouped by 
#'   haplotype.
#'   
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{privateAlleles}
#' 
#' @examples
#' data(msats.g)
#' 
#' propUniqueAlleles(msats.g)
#' 
#' @export
#' 
propUniqueAlleles <- function(g) { 
  sapply(locNames(g), function(locus) {
    id.a.freqs <- table(g@data$ids, g@data[[locus]])
    is.unique <- apply(id.a.freqs, 2, function(x) sum(x > 0) == 1)
    sum(is.unique)
  }) / numAlleles(g)
}