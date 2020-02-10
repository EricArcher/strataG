#' @title MAFFT Alignment
#' @description Align a set of sequences using the MAFFT executable.
#' 
#' @param x a list or a matrix of DNA sequences 
#'   (see \code{\link[ape]{write.dna}}).
#' @param run.label label for output alignment FASTA file.
#' @param delete.output logical. Delete output alignment FASTA file?
#' @param op gap opening penalty.
#' @param ep offset value, which works like gap extension penalty.
#' @param maxiterate number cycles of iterative refinement are performed.
#' @param quiet logical. Run MAFFT quietly?
#' @param num.cores The number of cores to use. If set to \code{NULL}, 
#'   the value will be what is reported by 
#'   \code{\link[parallel]{detectCores} - 1}. Passed to MAFFT argument 
#'   \code{--thread}.
#' @param opts character string other options to provide to command line.
#' @param simplify if \code{TRUE}, if \code{x} is a single sequence, 
#'   a single \code{DNAbin} object is returned, otherwise, a list of alignments
#'   is returned.
#' 
#' @note MAFFT is not included with \code{strataG} and must be downloaded 
#'   separately. Additionally, it must be installed such that it can be run from 
#'   the command line in the current working directory. See the vignette 
#'   for \code{external.programs} for installation instructions.
#' 
#' @return a \code{\link[ape]{DNAbin}} object with aligned sequences.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references 
#' Katoh, M., Kumar, M. 2002. MAFFT: a novel method for rapid multiple sequence 
#'   alignment based on fast Fourier transform. Nucleic Acids Res. 30:3059-3066.\cr 
#'   Available at: \url{http://mafft.cbrc.jp/alignment/software}
#' 
#' @examples \dontrun{
#' data(dolph.seqs)
#' dolph.aln <- mafft(dolph.seqs, op = 3, ep = 2)
#' dolph.aln
#' }
#' 
#' @export
#' 
mafft <- function(x, run.label = "align.mafft", delete.output = TRUE, 
                  op = 3, ep = 0.123, maxiterate = 0, quiet = FALSE, 
                  num.cores = 1, opts = "--auto", simplify = TRUE) {
  dna <- write.fasta(x, file = run.label) %>% 
    purrr::map(function(f) {
      aligned.fasta <- paste0("mafft.aligned", "_", f)
      mafft.cmd <- paste(
        "mafft", 
        opts, 
        "--op", op,
        "--ep", ep,
        "--maxiterate", maxiterate,
        ifelse(quiet, "--quiet", ""),
        "--thread", num.cores,
        f, ">", aligned.fasta
      )
      err.code <- system(mafft.cmd, intern = FALSE)
      if(!err.code == 0) return(NA)
      aligned <- read.fasta(aligned.fasta)
      file.remove(f)
      if(delete.output) file.remove(aligned.fasta)
      aligned
    })
  if(length(dna) == 1 & simplify) dna[[1]] else dna
}