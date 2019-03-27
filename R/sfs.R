#' @title Site Frequency Spectrum
#' @description Calculate the SFS from a data frame of SNP genotypes
#'   
#' @param x a data frame of SNP genotypes where the first two columns are id and
#'   strata designations and SNPs start on the third column. SNP genotypes are
#'   coded as integers where 0 and 2 are the major and minor homozygotes and 1
#'   is the heterozygote.
#' @param joint return the joint SFS if more than 1 strata is available?
#' 
#' @return Either a named vector of the 1D site frequency 
#'   spectrum if \code{joint = FALSE}, or a list of matrices of the 2D
#'   site frequency spectra for pairs of strata if \code{joint = TRUE}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
sfs <- function(x, joint = FALSE) {
  if(!joint | dplyr::n_distinct(x[, 2]) == 1) {
    count <- colSums(x[, -(1:2), drop = FALSE])
    loc.freq <- table(count)
    sfs <- stats::setNames(vector("numeric", nrow(x) + 1), 0:nrow(x))
    sfs[names(loc.freq)] <- loc.freq
    sfs
  } else {
    counts <- lapply(split(x[, -(1:2)], x[, 2]), colSums)
    n <- table(x[, 2])
    combn(unique(x[, 2]), 2, function(pair) {
      i <- pair[1]
      j <- pair[2]
      nrow <- (2 * n[i]) + 1
      ncol <- (2 * n[j]) + 1
      freqs <- table(counts[[i]], counts[[j]], dnn = pair) %>% 
        as.data.frame(stringsAsFactors = FALSE)
      sfs <- matrix(0, nrow = nrow, ncol = ncol)
      dimnames(sfs) <- stats::setNames(list(0:(nrow - 1), 0:(ncol - 1)), pair)
      for(r in 1:nrow(freqs)) sfs[freqs[r, 1], freqs[r, 2]] <- freqs[r, 3]
      t(sfs)
    }, simplify = FALSE)
  }
}