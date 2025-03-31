#' @title Convert FASTQ quality scores to base probability
#' @description Converts FASTQ quality scores encoded as ASCII characters to
#'   numeric probability that base was correctly called.
#' 
#' @param x a vector of quality scores represented as ASCII characters from 
#'   \code{!} to \code{~}. Each element in \code{x} can be a string of 
#'   characters or a single character.
#' @param phred.scale scale of PHRED encoding. Can be : 
#'   \code{"+33"} = range of \code{0:93} (e.g., Sanger, Illumina 1.8+), or
#'   \code{"+64"} = range of \code{-31:62} (e.g., Solexa, Illumina 1.3+, 1.5+)
#' @param simplify if there is only one set of quality scores, simplify 
#'   resulting list to a vector?
#' 
#' @return List of probabilities that each base represented by the characters in
#'   \code{x} was correctly called. List has one element for each set of 
#'   quality scores in \code{x}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples 
#' # each element is a single quality score
#' quality.1 <- sapply(sample(33:126, 15, rep = TRUE), intToUtf8)
#' prob.1 <- qual2prob(quality.1)
#' 
#' quality.1
#' prob.1
#' 
#' # each element is a string of quality scores
#' quality.2 <- replicate(5, {
#'   num.sites <- sample(5:10, 1)
#'   chars <- sapply(sample(33:126, num.sites, rep = TRUE), intToUtf8)
#'   paste(chars, collapse = "")
#' }, simplify = "vector")
#' prob.2 <- qual2prob(quality.2)
#' 
#' quality.2
#' prob.2 
#' 
#' @export
#' 
qual2prob <- function(x, phred.scale = c("+33", "+64"), simplify = TRUE) {
  if(!is.character(x)) stop("'x' must be a character vector")
  if(all(nchar(x) == 1)) x <- paste(x, collapse = "")
  x <- strsplit(x, "")
  
  phred.scale <- match.arg(phred.scale)
  result <- lapply(x, function(site) {
    sapply(site, function(char) {
      ascii <- utf8ToInt(char)
      if(!dplyr::between(ascii, 33, 126)) return(NA)
      q <- switch(phred.scale, '+33' = ascii - 33, '+64' = ascii - 64)
      1 - (10 ^ (-q / 10))
    }, USE.NAMES = FALSE)
  })
  
  if(simplify & length(result) == 1) result[[1]] else result
}