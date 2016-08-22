#' @title Quality Assurance/Quality Control
#' @description Conducts a suite of QA/QC tests. Summarizes missing data and 
#'   homozygosity by individual and locus, and looks for duplicate genotypes 
#'   (see \code{\link{dupGenotypes}}). For sequence data, identifies low 
#'   frequency substitutions (see \code{\link{lowFreqSubs}}), and computes 
#'   sequence likelihoods (see \code{\link{sequenceLikelihoods}}).
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param label optional label for output folder and prefix for files.
#' @param ... optional arguments to pass on to summary functions.
#' 
#' @return Files are written for by-sample and by-locus summaries, and duplicate 
#'   genotypes if any are found. If sequences are present, files are written 
#'   identifying low frequency substitutions and sequence likelihoods.\cr
#'   The return value is a list with the following elements:
#'   
#' \describe{
#'   \item{by.sample}{data.frame of by-sample summaries}
#'   \item{\code{by.locus}}{data.frame of by-locus summaries}
#'   \item{\code{dup.df}}{data.frame identifying potential duplicates}
#'   \item{\code{by.seq}}{list of low frequency substitutions and haplotype 
#'     likelihoods for each gene}
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{summarizeSamples}}, \code{\link{summarizeLoci}}, 
#'   \code{\link{dupGenotypes}}, \code{\link{lowFreqSubs}}, 
#'   \code{\link{sequenceLikelihoods}}
#' 
#' @importFrom utils write.csv
#' @export
#' 
qaqc <- function(g, label = NULL, ...) {
  cat("\n")
  cat(format(Sys.time()), ": Individual summaries\n")
  by.sample <- summarizeSamples(g)
  
  cat(format(Sys.time()), ": Locus summaries\n")
  by.locus <- summarizeLoci(g, TRUE)
  by.locus$All <- summarizeLoci(g)
  
  cat(format(Sys.time()), ": Duplicate genotypes\n")
  dup.df <- dupGenotypes(g, ...)
  
  # Sequence summaries
  by.seq <- if(!is.null(sequences(g))) {
    cat(format(Sys.time()), ": Sequence summaries\n")
    sapply(locNames(g), function(x) {
      x.seqs <- getSequences(sequences(g), loci = x, simplify = TRUE)
      list(
        low.freq.subs = lowFreqSubs(x.seqs, ...),
        seq.likelihoods = sequenceLikelihoods(x.seqs, ...)
      )
    }, simplify = FALSE)
  } else NULL
      
  cat(format(Sys.time()), ": Writing files\n")
  # Write summaries to files
  label <- if(is.null(label)) description(g) else label
  label <- gsub("[[:punct:]]", ".", label)
  fname <- paste(label, ".sample.summary.csv", sep = "")
  write.csv(by.sample, file = fname, row.names = FALSE)

  for(x in names(by.locus)) {
    fname <- paste(label, "locus.summary", x, "csv", sep = ".")
    by.locus[[x]] <- data.frame(
      locus = rownames(by.locus[[x]]), by.locus[[x]], stringsAsFactors = FALSE
    )
    write.csv(by.locus[[x]], file = fname, row.names = FALSE)
  }  
  
  if(!is.null(dup.df)) {
    fname <- paste(label, ".duplicate.samples.csv", sep = "")
    write.csv(dup.df, file = fname, row.names = FALSE)
  }
  
  if(!is.null(by.seq)) {
    for(x in names(by.seq)) {
      fname <- paste(label, "low.freq.subs", x, "csv", sep = ".")
      write.csv(by.seq[[x]]$low.freq.subs, file = fname, row.names = FALSE)
      fname <- paste(label, "sequence.likelihoods", x, "csv", sep = ".")
      write.csv(by.seq[[x]]$seq.likelihoods, file = fname, row.names = FALSE)
    }
  }
  
  cat("\n")
  invisible(list(by.sample = by.sample, by.locus = by.locus, dup.df = dup.df,
                 by.seq = by.seq))
}