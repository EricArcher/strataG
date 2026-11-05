#' @name iupac
#' @title IUPAC Codes
#' @description Calculate the correct IUPAC code(s) for a vector of nucleotides.
#'
#' @param bases character vector containing valid nucleotides or IUPAC codes.
#' @param all.codes logical. Return all valid codes (\code{TRUE}), or just the 
#'   most appropriate one (\code{FALSE}).
#' @param ignore.gaps logical. Ignore gaps at a site when computing the IUPAC codes.
#'
#' @return 
#' \tabular{ll}{
#'   \code{iupacCodes} \tab a character or vector of characters representing 
#'      the correct IUPAC code(s).\cr
#'   \code{iupacMat} \tab a logical matrix identifying valid IUPAC codes.\cr
#' }
#'
#' @author Eric Archer \email{eric.ivan.archer@@gmail.com}
#'
#' @examples
#' iupacCodes(c("a", "a", "g"))
#' 
#' iupacCodes(c("t", "c", "g"))
#'
#' iupacCodes(c("c", "t", "c", "c"), all = TRUE)
#' 
#' iupacCodes(c("c", "y", "c", "c"), all = TRUE)
#' 
#' iupacCodes(c("a", "g", "t", "a"), all = TRUE)
#' 
#' @export
#' 
iupacCodes <- function(bases, all.codes = FALSE, ignore.gaps = FALSE) {
  bases <- bases |> 
    as.character() |> 
    tolower()
  
  if(ignore.gaps) bases <- bases[!bases %in% c("-", ".")]
  
  base.rows <- which(rownames(iupac.mat) %in% bases)
  if(length(base.rows) == 0) stop("No valid IUPAC codes in 'bases'")
  
  valid.codes <- sapply(
    colnames(iupac.mat), 
    function(code) all(iupac.mat[base.rows, code])
  )
  
  codes <- colnames(iupac.mat)[valid.codes]
  if(all.codes) codes else codes[1]
}

#' @rdname iupac
#' @export
#' 
iupacMat <- function() iupac.mat
