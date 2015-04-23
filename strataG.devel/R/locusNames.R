#' @title Locus Names
#' @description Parse locus names from a matrix where alleles are in
#'   consecutive columns.
#'
#' @param loci a vector or matrix of alleles. If a matrix, each allele is in
#'   consecutive columns of the matrix.
#' @param ploidy integer specifying the ploidy
#' @param default.name name to provide for loci if no column names are provided
#'   in \code{loci}.
#'
#' @note If \code{loci} has column names, then columns of each locus must be
#'   named with the same roots. For example, diploid locus 'ABCD' would have
#'   two columns named something like 'ABCD.1' and 'ABCD.2', or
#'   'ABCD_A' and 'ABCD_B'. \cr
#'
#' @return A vector of locus names.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @export

locusNames <- function(loci, ploidy, default.name = "Locus") {
  if(is.vector(loci)) {
    loci <- cbind(loci)
    colnames(loci)
  }
  if(!is.null(colnames(loci))) return(colnames(loci))

  n.loc <- ncol(loci) / ploidy
  if(is.null(colnames(loci))) {
    # return generic names if no colnames assigned
    locus.nums <- formatC(1:n.loc, digits = floor(log10(n.loc)), flag = "0")
    paste(default.name, locus.nums, sep = "_")
  } else {
    locus.cols <- matrix(1:ncol(loci), nrow = ploidy)
    sapply(1:ncol(locus.cols), function(i) {
      loc.names <- colnames(loci)[locus.cols[, i]]
      max.length <- max(sapply(loc.names, nchar))
      # find location of first difference
      ptr <- 1
      while(length(unique(substr(loc.names, 1, ptr))) == 1 & ptr < max.length) ptr <- ptr + 1
      ptr <- ptr - 1
      if(ptr == 0) {
        loc.names[1] # return first column name if all different
      } else {
        # remove last same character if it is not alphanumeric
        if(!substr(loc.names[1], ptr, ptr) %in% c(LETTERS, letters, 0:9)) ptr <- ptr - 1
        substr(loc.names[1], 1, ptr)
      }
    })
  }
}