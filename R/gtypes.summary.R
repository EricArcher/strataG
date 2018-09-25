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
    dplyr::group_by(.data$stratum) %>% 
    dplyr::summarize(exptd.het = mean(.data$exptd.het, na.rm = TRUE))
  
  if(getPloidy(x) > 1) {
    het <- het %>% 
      dplyr::left_join(
        heterozygosity(x, TRUE, "obs") %>% 
          dplyr::group_by(.data$stratum) %>% 
          dplyr::summarize(obsvd.het = mean(.data$obsvd.het, na.rm = TRUE)),
        by = "stratum"
      )
  }
  
  strata.smry <- getNumInd(x, TRUE) %>% 
    dplyr::left_join(
      numMissing(x, TRUE) %>% 
        dplyr::group_by(.data$stratum) %>% 
        dplyr::summarize(num.missing = mean(.data$num.missing, na.rm = TRUE)),
      by = "stratum"
    ) %>% 
    dplyr::left_join(
      numAlleles(x, TRUE) %>% 
        dplyr::group_by(.data$stratum) %>% 
        dplyr::summarize(num.alleles = mean(.data$num.alleles, na.rm = TRUE)),
      by = "stratum"
    ) %>% 
    dplyr::left_join(
      propUniqueAlleles(x, TRUE) %>% 
        dplyr::group_by(.data$stratum) %>% 
        dplyr::summarize(prop.unique.alleles = mean(.data$prop.unique.alleles, na.rm = TRUE)),
      by = "stratum"
    ) %>% 
    dplyr::left_join(het, by = "stratum") %>% 
    as.data.frame()
  
  if(getPloidy(x) == 1) {
    strata.smry <- strata.smry %>% 
      dplyr::rename(
        num.haplotypes = .data$num.alleles,
        prop.unique.haplotypes = .data$prop.unique.alleles,
        haplotypic.diversity = .data$exptd.het
      )
  }
  
  x.seqs <- getSequences(x)
  seq.smry <- if(!is.null(x.seqs)) {
    do.call(rbind, lapply(names(x.seqs), function(gene) {
      x.seqs[[gene]] %>% 
        summarizeSeqs() %>% 
        as.data.frame() %>% 
        dplyr::summarize(
          num.seqs = dplyr::n(),
          mean.length = mean(.data$length),
          mean.num.ns = mean(.data$num.ns),
          mean.num.indels = mean(.data$num.indels)
        ) %>% 
        dplyr::mutate(locus = gene) %>% 
        dplyr::select(.data$locus, dplyr::everything()) %>% 
        as.data.frame() 
    }))
  } else NULL
        
  list(
    description = x@description,
    num.ind = getNumInd(x), 
    num.loc = getNumLoci(x), 
    num.strata = getNumStrata(x),
    schemes = if(!is.null(getSchemes(x))) colnames(getSchemes(x))[-1] else NULL,
    strata.smry = strata.smry,
    seq.smry = seq.smry
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
  if(!is.null(x$seq.smry)) {
    cat("\nSequence summary:\n")
    print(x$seq.smry)
  }
  cat("\n")
}