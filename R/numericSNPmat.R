#' @title Create a numeric SNP matrix
#' @description Create a matrix of SNPs coded as 0, 1, 2, where 0 and 2 are 
#'   the two homozygotes and 1 is the heterozygote.
#'
#' @param g a \linkS4class{gtypes} object. Must be diploid and have at least 
#'   one SNP with <= 2 alleles.
#' @param ref.allele an optional vector of reference alleles for each SNP. 
#'   If provided, it must be as long as there are biallelic SNPs in \code{g}. 
#'   If named, the SNP names must match those of all biallelic SNPs 
#'   in \code{g}. If set to \code{NULL} (default) the major allele at each 
#'   SNP is used as the reference.
#' 
#' @return a numeric matrix with all biallelic SNPs recoded to:  
#' 0 = homozygote for major or reference allele, 1 = heterozygote, 
#' 2 = homozygote for minor or alternate allele. Individual 
#' ids are in the rownames of the matrix.
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
numericSNPmat <- function(g, ref.allele = NULL) {
  if(getPloidy(g) != 2) stop("'g' must have diploid data")
  
  # find biallelic loci
  num.alleles <- numAlleles(g)
  to.keep <- num.alleles %>% 
    dplyr::filter(.data$num.alleles <= 2) %>% 
    dplyr::pull("locus")
  if(length(to.keep) == 0) stop("no biallelic loci found in 'g'")
  if(length(to.keep) != getNumLoci(g)) g <- g[, to.keep, ]
  
  # check reference allele vector
  ref.allele <- if(!is.null(ref.allele)) {
    if(length(ref.allele) != length(to.keep)) {
      stop("length of 'ref.allele' not equal to number of biallelic loci in 'g'")
    }
    # check for names
    if(!is.null(names(ref.allele))) {
      missing <- setdiff(to.keep, names(ref.allele))
      if(length(missing) != 0) {
        missing <- paste(missing, collapse = ", ")
        stop(
          "The following locus names in 'g' can't be found in 'ref.allele':",
          missing
        )
      }
      # filter provided vector for biallelic loci
      ref.allele[to.keep]
    } else {
      # add names to reference allele vector
      stats::setNames(ref.allele, to.keep)
    }
  } else {
    # use major alleles as reference
    sapply(alleleFreqs(g), function(x) names(which.max(x)))
  }

  # form matrix with one column per locus
  g.mat <- as.matrix(g, strata = FALSE, sep = "_<>_", one.col = TRUE)
  
  # count number of alleles for each individual that are not the reference
  mat <- sapply(colnames(g.mat[, -1, drop = FALSE]), function(loc) {
    loc.vec <- strsplit(g.mat[, loc, drop = FALSE], split = "_<>_")
    sapply(loc.vec, function(x) sum(x != ref.allele[loc]))
  })
  rownames(mat) <- g.mat[, 1]
  
  removed <- setdiff(getLociNames(g), to.keep)
  if(length(removed) > 0) {
    removed <- paste(removed, collapse = ", ")
    warning("The following loci have been removed: ", removed)
  }
  
  mat
}