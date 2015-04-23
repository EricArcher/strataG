#' @name summary,gtypes-method
#' @title Summarize gtypes Object
#' @description Generate a summary of a \code{gtypes} object.
#'  
#' @param object a \linkS4class{gtypes} object.
#' @param x list from summary.gtypes
#' @param ... other arguments (ignored).
#' 
#' @return a list with the following elements:
#' \tabular{ll}{
#'   \code{num.ind} \tab number of individuals.\cr
#'   \code{num.loc} \tab number of loci.\cr
#'   \code{num.strata} \tab number of strata.\cr
#'   \code{allele.freqs} \tab a list with tables of allele frequencies 
#'     by strata.\cr
#'   \code{strata.smry} \tab a by-strata data.frame summarizing haplotypes 
#'     or loci.\cr
#'   \code{seq.smry} \tab a summary of the sequence length and base 
#'     frequencies.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
setMethod("summary", "gtypes",
  function(object, ...) { 
    x <- object
  
    smry <- list(num.ind = nInd(x), num.loc = nLoc(x), num.strata = nStrata(x))
    smry$allele.freqs <- alleleFreqs(x, by.strata = TRUE)
    
    smry$strata.smry <- t(sapply(strataSplit(x), function(g) {
      if(g@ploidy == 1) {
        c(num.samples = nInd(g),
          num.missing = mean(numMissing(g)),
          num.haps = mean(numAlleles(g)), 
          hap.div = mean(obsvdHet(g)),
          pct.unique.haps = mean(pctUniqueHaps(g))
        )
      } else {
        c(num.samples = nInd(g),
          num.missing = mean(numMissing(g)),
          allelic.richness = mean(allelicRichness(g)),
          heterozygosity = mean(obsvdHet(g))
        )
      }
    }))
  
    smry$seq.smry <- if(!is.null(x@sequences)) {
      do.call(rbind, sapply(x@sequences@dna, function(dna) {
        dna.len <- unlist(lapply(dna, length))
        len.range <- range(dna.len)
        result <- data.frame(min.length = len.range[1], 
                             mean.length = round(mean(dna.len)), 
                             max.length = len.range[2]
        )
        cbind(result, rbind(base.freq(dna)))
      }, simplify = FALSE))
    } else NULL
  
    attr(smry, "description") <- x@description
    class(smry) <- c("gtypeSummary", "list")
    smry
})


#' @rdname summary-gtypes-method
#' @export
#' 
print.gtypeSummary <- function(x, ... ) { 
  ind.txt <- paste(x$num.ind, " sample", 
                   ifelse(x$num.ind > 1, "s", ""), sep = "")
  loc.txt <- paste(x$num.loc, " loc", 
                   ifelse(x$num.loc > 1, "i", "us"), sep = "")
  strata.txt <- paste(x$num.strata, " strat", 
                      ifelse(x$num.strata > 1, "a", "um"), sep = "")
  
  cat("\n")
  cat("<<<", attr(x, "description"), ">>>\n")
  cat("\nContents: ")
  cat(ind.txt, loc.txt, strata.txt, sep = ", ")
  cat("\n\nLoci: ")
  cat(names(x$allele.freqs), sep = ", ")
  cat("\n\nStrata summary:\n")
  print(x$strata.smry)
  if(!is.null(x$seq.smry)) {
    cat("\nSequence summary:\n")
    print(x$seq.smry)
  }
  cat("\n")
  
  invisible(x)
}