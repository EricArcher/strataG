#' @title Convert a data.frame to gtypes
#' @description Load allelic data from a data.frame or matrix into a 
#'   \linkS4class{gtypes} object. 
#' 
#' @param x a matrix or data.frame of genetic data.
#' @param ploidy number of number of columns in \code{x} storing alleles
#'   at each locus.
#' @param id.col column name or number where individual sample ids are stored.
#' @param strata.col column name or number where stratification is stored. If 
#'   NULL then all samples are in one (default) stratum.
#' @param loc.col column number of first allele of first locus.
#' @param sequences a list, matrix, \code{\link{DNAbin}}, or 
#'   \linkS4class{multidna} object containing sequences. 
#' @param description a label for the object (optional).
#' @param other a slot to carry other related information - unused in package
#'   analyses (optional).
#' 
#' @details
#' The genetic data in \code{x} starting at \code{loc.col} should be 
#' formatted such that every consecutive \code{ploidy} columns represent 
#' alleles of one locus. Locus names are taken from the column names in 
#' \code{x} and should be formatted with the same root locus name, with 
#' unique suffixes representing allels (e.g., for Locus1234: Locus1234.1 
#' and Locus1234.2, or Locus1234_A and Locus1234_B). \cr\cr
#' If sequences are provided in \code{sequences}, then they should be named 
#' and match haplotype labels in \code{loc.col} of \code{x}. If multiple 
#' genes are given as a \linkS4class{multidna}, then they should have the 
#' same names as column names in \code{x} from \code{loc.col} to the end.
#' 
#' @return a \linkS4class{gtypes} object.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' #--- create a diploid (microsatellite) gtypes object
#' data(dolph.msats)
#' data(dolph.strata)
#' msats.merge <- merge(dolph.strata[, c("ids", "fine")], dolph.msats, 
#'   all.y = TRUE)
#' msats.fine <- df2gtypes(msats.merge, ploidy = 2)
#' 
#' #' #--- create a haploid sequence (mtDNA) gtypes object
#' data(dolph.seqs)
#' seq.df <- dolph.strata[ c(1, 2, 1)]
#' colnames(seq.df)[3] <- "Haplotype"
#' dloop.broad <- df2gtypes(seq.df, ploidy = 1, sequences = dolph.seqs, 
#'   description = "dLoop: broad-scale stratification")
#' 
#' @export
#' 
df2gtypes <- function(x, ploidy, id.col = 1, strata.col = 2, loc.col = 3, 
                      sequences = NULL, description = NULL, other = NULL) {
  # check x
  if(!(is.matrix(x) | is.data.frame(x))) {
    stop("'x' must be a matrix or data.frame")
  }
  
  # extract locus and strata information
  rownames(x) <- x[, id.col]
  strata <- if(is.null(strata.col)) NULL else x[, strata.col]
  gen.data <- x[, loc.col:ncol(x), drop = FALSE]
  
  # return new gtypes object
  new("gtypes", gen.data = gen.data, ploidy = ploidy, strata = strata,
      sequences = sequences, description = description, other = other
  )
}