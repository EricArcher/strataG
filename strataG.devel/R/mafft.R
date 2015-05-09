#' @title MAFFT Alignment
#' 
#' @param x a \code{\link{gtypes}} object with aligned sequences or a list of aligned DNA sequences.
#' @param run.label label for output alignment FASTA file.
#' @param delete.output logical. Delete output alignment FASTA file?
#' @param op gap opening penalty.
#' @param ep offset value, which works like gap extension penalty.
#' @param maxiterate number cycles of iterative refinement are performed.
#' @param quiet logical. Run MAFFT quietly?
#' @param num.cores number of cores to be used. Passed to MAFFT argument \code{--thread}.
#' @param opts character string other options to provide to command line.
#' 
#' @note Formats and executes a call to the executable \code{mafft}, assuming that is installed
#'   on the system and available at the command line. 
#' 
#' @return list of aligned sequences
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references MAFFT is available at: \url{http://mafft.cbrc.jp/alignment/software}
#' 
#' @export
#' 
mafft <- function(x, run.label = "align.mafft", delete.output = TRUE, 
                        op = 3, ep = 0.123, maxiterate = 0, quiet = FALSE, 
                        num.cores = 1, opts = "--auto") {
  
  in.fasta <- paste(run.label, ".in.fasta", sep = "")
  aligned.fasta <- paste(run.label, ".aligned.fasta", sep = "")
  write.fasta(x, file = in.fasta)
  mafft.cmd <- paste("mafft", 
                     opts, 
                     "--op", op,
                     "--ep", ep,
                     "--maxiterate", maxiterate,
                     ifelse(quiet, "--quiet", ""),
                     "--thread", num.cores,
                     in.fasta, ">", aligned.fasta
  )
  err.code <- system(mafft.cmd, intern = FALSE)
  if(!err.code == 0) return(NA)
  aligned <- read.fasta(aligned.fasta)
  file.remove(in.fasta)
  if(delete.output) file.remove(aligned.fasta)
  aligned
}