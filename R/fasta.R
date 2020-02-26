#' @name fasta
#' @aliases fasta
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
#' @export
#' 
read.fasta <- function(file) {
  ape::read.dna(
    file, format = "fasta", as.character = TRUE, as.matrix = FALSE
  ) %>% 
    # replace ?s with Ns and convert to lower-case
    purrr::map(function(x) tolower(gsub("\\?", "n", x))) %>% 
    ape::as.DNAbin() %>% 
    as.list()
}


#' @rdname fasta
#' @export
#' 
write.fasta <- function(x, file = NULL) {
  fasta.func <- function(dna, f) {
    ape::write.dna(
      dna, file = f, format = "fasta", nbcol = -1, 
      colsep = "", indent = 0, blocksep = 0
    )
    f
  }
  
  folder <- if(is.null(file)) "." else dirname(file)
  if(is.gtypes(x)) {
    if(!is.null(file)) file <- basename(file)
    file <- .getFileLabel(x, file)
  }
  if(is.null(file)) file <- "sequences.fasta"
  if(!grepl("((\\.fas)$)|((\\.fasta)$)", tolower(file))) {
    file <- paste0(file, ".fasta")
  }
  
  fname <- if(inherits(x, "multidna") | is.gtypes(x)) {
    x <- apex::getSequences(as.multidna(x), simplify = TRUE)
    sapply(names(x), function(gene) {
      fasta.func(x[[gene]], file.path(folder, paste0(gene, " ", file)))
    })
  } else {
    x <- if(inherits(x, "DNAbin")) as.character(x) else ape::as.DNAbin(x)
    fasta.func(x, file.path(folder, file))
  }
  invisible(fname)
}
