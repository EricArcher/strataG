#' @title Trim N's From Sequences
#' @description Removes N's from beginning and end of sequences.
#' 
#' @param x a list of DNA sequences.
#' 
#' @return a list of DNA sequences.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
trimNs <- function(x) {
  x <- as.multidna(x)
  result <- lapply(x@dna, function(dna) {
    dna <- as.character(dna)
    lapply(1:nrow(dna), function(i) {
      seq.vec <- paste(dna[i, ], collapse = "")
      start <- gregexpr("^[n]+", seq.vec)[[1]]
      end <- gregexpr("[n]+$", seq.vec)[[1]]
      start <- ifelse(start == -1, 1, attr(start, "match.length") + 1)
      end <- ifelse(end == -1, nchar(seq.vec), end)
      dna[i, start:end]
    })
  })
  
  if(length(result) == 1) {
    as.DNAbin(result[[1]])
  } else {
    names(result) <- names(x@dna)
    new("multidna", result)
  }
}