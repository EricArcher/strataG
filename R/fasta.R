#' @title Read and Write FASTA
#' @description Read and write FASTA formatted files of sequences.
#' 
#' @param file a FASTA-formatted file of sequences.
#' @param x a list or a matrix of DNA sequences (see \code{\link[ape]{write.dna}}), 
#'   or a \code{\link{gtypes}} object with sequences.
#' 
#' @return \describe{
#'   \item{read.fasta}{a set of sequences in DNAbin format}
#'   \item{write.fasta}{invisbly, name(s) of file(s) written}
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @name fasta
#' @aliases fasta
#' @export
#' 
read.fasta <- function(file) {
  dna.seq <- read.dna(file, format = "fasta", as.character = TRUE, 
                      as.matrix = FALSE)
  # replace ?s with Ns and convert to lower-case
  as.list(as.DNAbin(lapply(dna.seq, function(x) tolower(gsub("\\?", "n", x)))))
}


#' @rdname fasta
#' @export
#' 
write.fasta <- function(x, file = "sequences.fasta") {
  if(is.gtypes(x)) x <- sequences(x)
  fname <- if(inherits(x, "multidna")) {
    x <- getSequences(x, simplify = FALSE)
    sapply(names(x), function(gene) {
      write.dna(x[[gene]], file = paste(gene, file), format = "fasta", 
                nbcol = -1, colsep = "", indent = 0, blocksep = 0)
    })
  } else {
    if(inherits(x, "DNAbin")) x <- as.character(x)
    write.dna(x, file = file, format = "fasta", nbcol = -1, colsep = "", 
              indent = 0, blocksep = 0)
    file
  }
  invisible(file)
}
