#' @title Trim N's From Sequences
#' @description Removes N's from beginning and end of sequences.
#' 
#' @param x a \code{\link[ape]{DNAbin}} object or list or matrix that can be
#'   coerced into one.
#' 
#' @return sequences with beginning and trailing N's removed.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#'  test.seqs <- list(
#'    A = c(rep("n", 5), "a", "c", "g", "t", rep("n", 3)),
#'    B = c(rep("n", 3), "a", "c", "g", "t", rep("n", 5)),
#'    C = c("a", "c", "g", "t", rep("n", 8))
#'  )
#' 
#' test.seqs
#' trimmed <- trimNs(test.seqs)  
#' as.character(trimmed)
#'  
#' @export
#' 
trimNs <- function(x) {
  if(!inherits(x, "DNAbin")) {
    if(inherits(x, "list") | inherits(x, "matrix")) {
      x <- ape::as.DNAbin(x)
    } else {
      stop("'x' must be a DNAbin object or a list or matrix of sequences")
    }
  }
  
  dna <- as.character(as.list(x))
  result <- lapply(1:length(dna), function(i) {
    seq.vec <- paste(dna[[i]], collapse = "")
    start <- gregexpr("^[n]+", seq.vec)[[1]]
    end <- gregexpr("[n]+$", seq.vec)[[1]]
    start <- ifelse(start == -1, 1, attr(start, "match.length") + 1)
    end <- ifelse(end == -1, nchar(seq.vec), end - 1)
    dna[[i]][start:end]
  })
  
  stats::setNames(ape::as.DNAbin(result), names(x))
}