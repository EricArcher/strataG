#' @title Show a gtypes object
#' @description Show a gtypes object
#' 
#' @param object a \linkS4class{gtypes} object.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
methods::setMethod("show", "gtypes", function(object) {
  .printBaseSmry(.baseSmry(object))
  invisible(NULL)
})


#' @keywords internal
#' 
.baseSmry <- function(object) {
  # object is a gtypes object
  strata.smry <- getNumInd(object, TRUE) |> 
    dplyr::left_join(
      numMissing(object, TRUE) |> 
        dplyr::group_by(.data$stratum) |> 
        dplyr::summarize(num.missing = mean(.data$num.missing, na.rm = TRUE)),
      by = "stratum"
    ) |> 
    dplyr::left_join(
      numAlleles(object, TRUE) |> 
        dplyr::group_by(.data$stratum) |> 
        dplyr::summarize(num.alleles = mean(.data$num.alleles, na.rm = TRUE)),
      by = "stratum"
    ) |> 
    as.data.frame()
  
  if(getPloidy(object) == 1) {
    strata.smry <- strata.smry |> 
      dplyr::rename(num.haplotypes = .data$num.alleles)
  }
  
  g.seqs <- getSequences(object, as.haplotypes = TRUE, simplify = FALSE)
  seq.smry <- if(!is.null(g.seqs)) {
    do.call(rbind, lapply(names(g.seqs), function(gene) {
      g.seqs[[gene]] |> 
        summarizeSeqs() |> 
        as.data.frame() |> 
        dplyr::summarize(
          num.seqs = dplyr::n(),
          mean.length = mean(.data$length),
          mean.num.ns = mean(.data$num.ns),
          mean.num.indels = mean(.data$num.indels)
        ) |> 
        dplyr::mutate(locus = gene) |> 
        dplyr::select(.data$locus, dplyr::everything()) |> 
        as.data.frame() 
    }))
  } else NULL
  
  list(
    description = getDescription(object),
    num.ind = getNumInd(object), 
    num.loc = getNumLoci(object), 
    num.strata = getNumStrata(object),
    schemes = if(!is.null(getSchemes(object))) {
      colnames(getSchemes(object))[-1] 
    } else NULL,
    strata.smry = strata.smry,
    seq.smry = seq.smry, 
    other = names(getOther(object))
  )
}


#' @keywords internal
#' 
.printSmryHeader <- function(x) {
  # x is the output from .baseSmry()
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
}


#' @keywords internal
#' 
.printBaseSmry <- function(x) {
  # x is the output from .baseSmry()
  .printSmryHeader(x)
  if(!is.null(x$schemes)) {
    cat("\nStratification schemes:", paste(x$schemes, collapse = ", "))
  }
  if(!is.null(x$other)) {
    cat("\nOther info:", paste(x$other, collapse = ", "))
  }
  cat("\n\nStrata summary:\n")
  print(x$strata.smry)
  if(!is.null(x$seq.smry)) {
    cat("\nSequence summary:\n")
    print(x$seq.smry)
  }
  cat("\n")
}
