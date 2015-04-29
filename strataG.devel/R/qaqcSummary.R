#' @name qaqcSummary
#' 
#' @title QA/QC Summaries
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param label optional label for output folder and prefix for files.
#' @param ... optional arguments to pass on to summary functions.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
qaqcSummary <- function(g, label = NULL, ...) {
  cat("\n")
  cat(format(Sys.time()), ": Individual summaries\n")
  num.loc <- nLoc(g)
  by.sample <- do.call(rbind, lapply(indNames(g), function(id) {
    smry <- sapply(locNames(g), function(loc) {
      gt <- genotype(id, loc, g)
      missing <- any(is.na(gt))
      hmzgs <- if(missing) NA else length(unique(gt)) == 1
      c(missing = missing, hmzgs = hmzgs)
    })
    
    missing <- sum(smry["missing", ], na.rm = TRUE)
    c(sample = id, 
      num.loci.missing.genotypes = missing,
      pct.loci.missing.genotypes = missing / num.loc,
      pct.loci.homozygous = mean(smry["hmzgs", ], na.rm = TRUE)
    )
  }))
  
  cat(format(Sys.time()), ": Locus summaries\n")
  by.locus <- summarizeLoci(g, TRUE)
  by.locus$All <- summarizeLoci(g)
  
  cat(format(Sys.time()), ": Duplicate genotypes\n")
  dup.df <- dupGenotypes(g, ...)
  
  # Sequence summaries
  by.seq <- if(!is.null(sequences(g))) {
    cat(format(Sys.time()), ": Sequence summaries\n")
    sapply(seqNames(g), function(x) {
      list(
        low.freq.subs = lowFreqSubs(sequences(g, x), ...),
        hap.likelihoods = haplotypeLikelihoods(sequences(g, x), ...)
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
      fname <- paste(label, "haplotype.likelihoods", x, "csv", sep = ".")
      hl <- cbind(delta.logLik = by.seq[[x]]$hap.likelihoods)
      write.csv(hl, file = fname)
    }
  }
  
  cat("\n")
  invisible(list(by.sample = by.sample, by.locus = by.locus, dup.df = dup.df,
                 by.seq = by.seq))
}