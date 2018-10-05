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
#' @examples 
#' data(bowhead.snps)
#' snps <- df2gtypes(bowhead.snps, ploidy = 2)
#' 
#' num.mat <- numericSNPmat(snps)
#' head(num.mat[, 1:5])
#'
#' @export
#' 
numericSNPmat <- function(g) {
  if(getPloidy(g) != 2) stop("'g' must have diploid data")
  
  nums <- c('1.1' = 0, '1.2' = 1, '2.2' = 2)
  
  g.mat <- as.matrix(g, sep = "_", one.col = TRUE)
  
  mat <- apply(g.mat[, -(1:2)], 2, function(loc) {
    loc <- strsplit(loc, split = "_")
    lvls <- sort(unique(unlist(loc)))
    lvls <- lvls[!is.na(lvls)]
    if(length(lvls) > 2) return(rep(NA, nrow(g.mat)))
    unname(sapply(loc, function(x) nums[paste(match(x, lvls), collapse = ".")]))
  })
  rownames(mat) <- g.mat[, 1]
  to.delete <- which(apply(mat, 2, function(x) all(is.na(x))))
  if(length(to.delete) > 0) {
    mat <- mat[, -to.delete]
    to.delete <- paste(colnames(mat)[to.delete], collapse = ", ")
    warning("The following loci have been removed: ", to.delete)
  }
  mat
}