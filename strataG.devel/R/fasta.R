#' @name fasta
#' @title Read and Write FASTA
#' @description Read and write FASTA formatted files of sequences.
#' 
#' @param file a FASTA-formatted file of sequences.
#' @param x a list of DNA sequences or a haploid \code{\link{gtypes}} object 
#'   with sequences. 
#' 
#' @return from \code{read.fasta}, a list of DNA sequences.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @aliases fasta
#' @export
#' 
read.fasta <- function(file) {
  dna.seq <- read.dna(file, format = "fasta", as.character = TRUE, 
                      as.matrix = FALSE)
  # replace ?s with Ns and convert to lower-case
  as.DNAbin(lapply(dna.seq, function(x) tolower(gsub("\\?", "n", x))))
}

#' @rdname fasta
#' @export
#' 
write.fasta <- function(x, file = "sequences.fasta") {
  write.dna(x, file = file, format = "fasta", nbcol = -1, colsep = "", 
            indent = 0, blocksep = 0)
}
