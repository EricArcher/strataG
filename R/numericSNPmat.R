#' @title Create a numeric SNP matrix
#' @description Create a matrix of SNPs coded as 0, 1, 2, where 0 and 2 are 
#'   the two homozygotes and 1 is the heterozygote.
#'
#' @param g a \linkS4class{gtypes} object.
#' 
#' @return a data.frame with all biallelic loci recoded.
#' 
#' @note \code{g} must be a diploid \code{gtypes} object and have at least one 
#'   locus with one or two alleles.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export

numericSNPmat <- function(g) {
  if(ploidy(g) != 2) stop("'g' must have diploid data")
  
  num.alleles <- numAlleles(g)
  biallelic <- names(num.alleles)[num.alleles <= 2]
  if(length(biallelic) == 0) {
    warning("No loci are biallelic. No file written.")
    return(NULL)
  }
  
  g <- g[, biallelic, ]
  mat <- sapply(locNames(g), function(locus) {
    arr <- as.array(g, locus = locus, drop = TRUE)
    apply(arr, 1, function(x) {
      switch(paste(sort(x), collapse = "."), '1.1' = 0, '1.2' = 1, '2.2' = 2, NA)
    })
  })
  rownames(mat) <- indNames(g)
  mat
}