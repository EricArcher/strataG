#' @title Summarize gtypes Object
#' @description Generate a summary of a \code{gtypes} object.
#'  
#' @param object a \linkS4class{gtypes} object.
#' @param x list from summary.gtypes
#' @param ... other arguments (ignored).
#' 
#' @return a list with the following elements:
#' \describe{
#'   \item{\code{num.ind}}{number of individuals}
#'   \item{\code{num.loc}}{number of loci}
#'   \item{\code{num.strata}}{number of strata}
#'   \item{\code{unstratified}}{number of unstratified samples}
#'   \item{\code{schemes}}{names of stratification schemes}
#'   \item{\code{allele.freqs}}{a list with tables of allele frequencies by strata}
#'   \item{\code{strata.smry}}{a by-strata data.frame summarizing haplotypes or loci}
#'   \item{\code{locus.smry}}{a data.frame summarizing each locus for 
#'     non-haploid objects, \code{NULL} for haploid objects}
#'   \item{\code{seq.smry}}{a summary of the sequence length and base frequencies}
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @name summary,gtypes-method
#' @aliases summary summary.gtypes
#' 
NULL

.baseSmry <- function(x) {
  het <- heterozygosity(x, TRUE, "exp") %>% 
    dplyr::group_by(stratum) %>% 
    dplyr::summarize(exptd.het = mean(exptd.het, na.rm = TRUE))
  
  if(ploidy(x) > 1) {
    het <- het %>% 
      left_join(
        heterozygosity(x, TRUE, "obs") %>% 
          dplyr::group_by(stratum) %>% 
          dplyr::summarize(obsvd.het = mean(obsvd.het, na.rm = TRUE)),
        by = "stratum"
      )
  }
  
  strata.smry <- getNumInd(x, TRUE) %>% 
    left_join(
      numMissing(x, TRUE) %>% 
        dplyr::group_by(stratum) %>% 
        dplyr::summarize(num.missing = mean(num.missing, na.rm = TRUE)),
      by = "stratum"
    ) %>% 
    left_join(
      numAlleles(x, TRUE) %>% 
        dplyr::group_by(stratum) %>% 
        dplyr::summarize(num.alleles = mean(num.alleles, na.rm = TRUE)),
      by = "stratum"
    ) %>% 
    left_join(
      propUniqueAlleles(x, TRUE) %>% 
        dplyr::group_by(stratum) %>% 
        dplyr::summarize(prop.unique.alleles = mean(prop.unique.alleles, na.rm = TRUE)),
      by = "stratum"
    ) %>% 
    left_join(het, by = "stratum") %>% 
    as.data.frame()
  
  if(ploidy(x) == 1) {
    strata.smry <- strata.smry %>% 
      rename(
        num.haplotypes = num.alleles,
        prop.unique.haplotypes = prop.unique.alleles,
        haplotypic.diversity = exptd.het
      )
  }
        
  list(
    description = x@description,
    num.ind = getNumInd(x), 
    num.loc = getNumLoci(x), 
    num.strata = getNumStrata(x),
    schemes = if(!is.null(schemes(x))) colnames(schemes(x))[-1] else NULL,
    strata.smry = strata.smry
  )
}

.printBaseSmry <- function(x) {
  ind.txt <- paste(x$num.ind, " sample", 
                   ifelse(x$num.ind > 1, "s", ""), sep = "")
  loc.txt <- paste(x$num.loc, " loc", 
                   ifelse(x$num.loc > 1, "i", "us"), sep = "")
  strata.txt <- paste(x$num.strata, " strat", 
                      ifelse(x$num.strata > 1, "a", "um"), sep = "")
  
  cat("\n")
  cat("<<<", x$description, ">>>\n")
  cat("\nContents: ")
  cat(ind.txt, loc.txt, strata.txt, sep = ", ")
  if(!is.null(x$schemes)) cat("\nStratification schemes:", paste(x$schemes, collapse = ", "))
  cat("\n\nStrata summary:\n")
  print(x$strata.smry)
  cat("\n")
}


#' #' @rdname summary-gtypes-method
#' #' @export
#' #' 
#' setMethod("summary", "gtypes", function(object, type = c("loci", "individual"), ...) { 
#'   smry <- .baseSmry(object)
#'   smry$allele.freqs <- alleleFreqs(object, by.strata = TRUE)
#'   smry$sample.smry <- summarizeSamples(object)
#'   smry$locus.smry <- if(ploidy(object) > 1) summarizeLoci(object) else NULL
#'   smry$seq.smry <- if(!is.null(sequences(object))) {
#'     sequences <- getSequences(sequences(object), simplify = FALSE)
#'     do.call(rbind, sapply(sequences, function(dna) {
#'       dna <- as.matrix(dna)
#'       dna.len <- unlist(lapply(dna, length))
#'       len.range <- range(dna.len)
#'       result <- data.frame(
#'         num.seqs = nrow(dna), 
#'         min.length = len.range[1], 
#'         mean.length = round(mean(dna.len)), max.length = len.range[2]
#'       )
#'       cbind(result, rbind(base.freq(dna)))
#'     }, simplify = FALSE))
#'   } else NULL
#'   
#'   class(smry) <- c("gtypeSummary", "list")
#'   smry
#' })
#' 
#' 
#' #' @rdname summary-gtypes-method
#' #' @export
#' #' 
#' print.gtypeSummary <- function(x, ... ) { 
#'   .printBaseSmry(x)
#'   if(!is.null(x$locus.smry)) {
#'     cols <- c(1, 3, 5, 7)
#'     num.rows <- nrow(x$locus.smry)
#'     if(num.rows > 20) {
#'       cat("\nLocus summary (first and last 10):\n")
#'       print(x$locus.smry[1:10, cols, drop = FALSE])
#'       cat("---\n")
#'       print(x$locus.smry[(num.rows - 10):num.rows, cols, drop = FALSE])
#'     } else {
#'       cat("\nLocus summary:\n")
#'       print(x$locus.smry[, cols, drop = FALSE])
#'     }
#'   }
#'   if(!is.null(x$seq.smry)) {
#'     cat("\nSequence summary:\n")
#'     print(x$seq.smry)
#'   }
#'   cat("\n")
#'   invisible(x)
#' }