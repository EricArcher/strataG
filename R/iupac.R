#' @name iupac
#' @title IUPAC Codes
#' @description Calculate the correct IUPAC code for a vector of nucleotides.
#'
#' @param bases character vector containing valid nucleotides or IUPAC codes.
#' @param ignore.gaps logical. Ignore gaps at a site when creating consensus. 
#'   If true, then bases with a gap are removed before consensus is calculated. 
#'   If false and a gap is present, then the result is a gap.
#'
#' @return 
#' \tabular{ll}{
#'   \code{iupacCode} \tab a character representing the correct IUPAC code 
#'     \code{bases}.\cr
#'   \code{validIupacCodes} \tab a character vector of all valid IUPAC 
#'     codes for \code{bases}.\cr
#'   \code{iupacMat} \tab a logical matrix identifying valid IUPAC codes.\cr
#' }
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @seealso \code{\link{validIupacCodes}}
#' 
#' @examples
#' iupacCode(c("a", "a", "g"))
#' 
#' iupacCode(c("t", "c", "g"))
#'
#' validIupacCodes(c("c", "t", "c", "c"))
#' 
#' validIupacCodes(c("c", "y", "c", "c"))
#' 
#' validIupacCodes(c("a", "g", "t", "a"))
#' 
#' @export
#' 
iupacCode <- function(bases, ignore.gaps = FALSE) {
  bases <- as.character(bases)
  if(ignore.gaps) bases <- bases[!bases %in% c("-", ".")]
  validIupacCodes(bases)[1]
}

#' @rdname iupac
#' @export
#' 
validIupacCodes <- function(bases) {
  bases <- tolower(bases)
  base.rows <- which(rownames(iupac.mat) %in% bases)
  if(length(base.rows) == 0) stop("No valid IUPAC codes in 'bases'")
  valid.codes <- sapply(
    colnames(iupac.mat), 
    function(code) all(iupac.mat[base.rows, code])
  )
  colnames(iupac.mat)[valid.codes]
}

#' @rdname iupac
#' @export
#' 
iupacMat <- function() iupac.mat