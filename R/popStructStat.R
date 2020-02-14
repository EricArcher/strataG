#' @title Population structure statistics
#' @description Population structure statistics
#'
#' @param g a \linkS4class{gtypes} object.
#' @param nrep number specifying number of permutation replicates to use for 
#'   permutation test.
#' @param keep.null logical. Keep the null distribution from the 
#'   permutation test?
#' @param prime.type type of G'st to calculate. Can be "hedrick" or "nei".
#' @param model,gamma,pairwise.deletion parameters passed to 
#'   \code{\link[ape]{dist.dna}}. Note that defaults for these arguments 
#'   (in particular \code{model}) are the same as in \code{dist.dna}.
#' @param ... optional arguments passed to or from other functions.
#'
#' @return A list with three elements: \describe{
#'   \item{stat.name}{the name of the statistic.}
#'   \item{result}{a vector of the statistic estimate and the p-value, 
#'     if replicates were conducted.}
#'   \item{null.dist}{a vector of the null distribution from the permutations.}
#' }
#' 
#' @note If \code{strata.mat} is provided, it must be a numeric matrix of 
#'   integers from \code{0} to \code{k - 1}, where \code{k} is the number of 
#'   strata. Each column is a separate permutation and the first column is 
#'   assumed to represent the original stratification. If not provided 
#'   (\code{strata.mat = NULL}), stratification is taken from \code{g}. 
#'   This argument is primarily used internally by \code{\link{popStructTest}}.
#'
#' @references \describe{
#'   \item{Hstats, statFis}{Nei, M. and R.K. Chesser. 1983. Estimation of
#'     fixation indices and gene diversities. Ann. Hum. Genet. 47:253-259.}
#'   \item{statFst}{Weir, B. and Cockerham, C. 1984. Estimating F-Statistics 
#'     for the Analysis of Population Structure. Evolution 38(6):1358-1370. 
#'     doi:10.2307/2408641}
#'   \item{statGst, statGstPrime}{Nei, M. 1973. Analysis of gene diversity in 
#'     subdvidided populations. Proc. Nat. Acad. Sci. 70(12):3321-3323.  
#'     Hedrick, P.W. 2005. A standardized genetic differentiation measure. 
#'     Evolution 59(8):1633-1638.}
#'   \item{statFstPrime, statGstDblPrime}{Miermans, P.G. and P.W. Hedrick. 2011. 
#'     Assessing population structure: FST and related measures. Molecular 
#'     Ecology Resources 11:5-18.}
#'   \item{statJostD}{Jost, L. 2008. GST and its relatives do not measure 
#'     differentiation. Molecular Ecology 17:4015-4026.}
#'   \item{statPHIst}{Excoffier, L., P.E. Smouse, and J.M. Quattro. 1992. 
#'     Analysis of molecular variance inferred from metric distances among 
#'     DNA haplotypes: Application to human mitochondrial DNA restriction data. 
#'     Genetics 131:479-491.}
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @name popStructStat
#' @useDynLib strataG, .registration = TRUE
#'
NULL

#----- Supporting Functions -----

#' @noRd
#' 
# create a list of formatted numeric loci and strata matrices for input to C
#   population structure functions
.formatCinput <- function(g, nrep, keep.null, hap.dist = FALSE, ...) {
  # remove unstratified samples
  g <- g[, , getStrataNames(g)]
  
  # delete loci with no genotypes in at least one stratum
  freqs <- alleleFreqs(g, TRUE)
  to.delete <- sapply(freqs, function(loc) {
    !all(apply(loc, 2, function(st) any(loc > 0)))
  })
  if(sum(to.delete) > 0) {
    warning(paste(
      "The following ", sum(to.delete), 
      " loci will be removed because they have no genotypes in one or more strata: ",
      paste(names(freqs)[to.delete], collapse = ", ")
    ))
    if(sum(!to.delete) == 0) stop("no loci available for analysis.")
    g <- g[, names(freqs)[!to.delete], ]
  }
  if(getNumStrata(g) < 2) stop("'g' must have more than one stratum defined.")
  
  # create matrix of numeric loci
  input <- list(
    loci =.stackedAlleles(g, alleles2integer = TRUE) %>% 
      dplyr::select(-.data$stratum, -.data$allele)
  )
  
  # create matrix of permuted numeric strata
  st <- getStrata(g)
  if(any(is.na(st))) stop("cannot run with unstratified samples")
  strata.num <- cbind(as.numeric(factor(st)) - 1)
  rownames(strata.num) <- names(st)
  if(is.null(nrep)) nrep <- 0
  strata <- if(nrep < 1) {
    strata.num
  } else {
    cbind(strata.num, sapply(1:nrep, function(i) sample(strata.num)))
  }
  input$strata <- strata[unique(input$loci$id), , drop = FALSE]
  
  input$ploidy <- getPloidy(g)
  input$keep.null = keep.null
  
  func.args <- lapply(
    as.list(match.call()[-1]), 
    eval, 
    envir = parent.frame()
  )
  
  input$prime.type <- if(is.null(func.args$prime.type)) {
    1
  } else {
    switch(
      match.arg(func.args$prime.type, choices = c("nei", "hedrick")), 
      nei = 0, 
      hedrick = 1
    )
  }
  
  if(hap.dist) {
    if(is.null(func.args$model)) func.args$model <- "K80"
    if(is.null(func.args$gamma)) func.args$gamma <- FALSE
    if(is.null(func.args$pairwise.deletion)) func.args$pairwise.deletion <- TRUE
    alleles <- getAlleleNames(g)
    make.unit.dist <- is.null(getSequences(g))
    input$hap.dist <- sapply(colnames(input$loci)[-1], function(gene) {
      haps <- alleles[[gene]]
      if(make.unit.dist) {
        # format distances for Fst (all 1s, and 0s on the diagonal)
        hd <- matrix(1, nrow = length(haps), ncol = length(haps), 
                     dimnames = list(haps, haps))
        diag(hd) <- 0
        hd
      } else {
        ape::dist.dna(
          getSequences(g)[[gene]][haps], 
          model = func.args$model, 
          gamma = func.args$gamma, 
          pairwise.deletion = func.args$pairwise.deletion, 
          as.matrix = TRUE
        )[haps, haps]
      }
    }, USE.NAMES = TRUE, simplify = FALSE)
  }
  
  input$loci <- input$loci %>% 
    dplyr::select(-.data$id) %>% 
    as.matrix
  input
}


#' @noRd
#' 
# create returned list from result vector
.formatResult <- function(stat.name, result = NULL, keep.null = FALSE) {
  if(is.null(result)) {
    return(list(
      stat.name = stat.name,
      result = c(estimate = NA, p.val = NA),
      null.dist = NULL
    ))
  }
  if(length(result) == 1) keep.null <- FALSE
  p.val <- if(length(result) == 1) NA else mean(result >= result[1], na.rm = TRUE) 
  if(is.nan(p.val)) p.val <- NA
  return(list(
    stat.name = stat.name, 
    result = c(estimate = result[1], p.val = p.val),
    null.dist = if(keep.null) result[-1] else NULL
  ))
}


#----- Population Structure Functions -----

#' @noRd
#' 
.statChi2 <- function(input) {
  if(is.null(input)) return(.formatResult("Chi2", NULL, input$keep.null))
  result <- statChi2_C(input$loci, input$strata)
  .formatResult("Chi2", result, input$keep.null)
}

#' @rdname popStructStat
#' @export
#' 
statChi2 <- function(g, nrep = NULL, keep.null = FALSE, ...) {
  .statChi2(.formatCinput(g, nrep, keep.null))
}

#-----


#' @noRd
#' 
.statJostD <- function(input) {
  if(is.null(input)) return(.formatResult("D", NULL, input$keep.null))
  if(input$ploidy == 1) return(.formatResult("D", NULL, input$keep.null))
  result <- statJostD_C(input$loci, input$strata)
  .formatResult("D", result, input$keep.null)
}

#' @rdname popStructStat
#' @export
#' 
statJostD <- function(g, nrep = NULL, keep.null = FALSE, ...) {
  .statJostD(.formatCinput(g, nrep, keep.null))
}

#-----


#' @noRd
#' 
.statFis <- function(input) {
  if(is.null(input)) return(.formatResult("Fis", NULL, input$keep.null))
  if(input$ploidy == 1) return(.formatResult("Fis", NULL, input$keep.null))
  result <- statFis_C(input$loci, input$strata)
  .formatResult("Fis", result, input$keep.null)
}

#' @rdname popStructStat
#' @export
#' 
statFis <- function(g, nrep = NULL, keep.null = FALSE, ...) {
  .statFis(.formatCinput(g, nrep, keep.null))
}

#-----


#' @noRd
#' 
.statFst <- function(input) {
  if(is.null(input)) return(.formatResult("Fst", NULL, input$keep.null))
  # for haploid data, run Phist with unity haplotype distance matrix
  result <- if(input$ploidy == 1) {
    input$hap.dist <- lapply(1:ncol(input$loci), function(i) {
      num.haps <- max(input$loci[, i], na.rm = TRUE) + 1
      hd <- matrix(1, nrow = num.haps, ncol = num.haps)
      diag(hd) <- 0
      hd
    })
    result <- .statPhist(input)
    result$stat.name <- "Fst"
    return(result)
  }
  result <- statFst_C(input$loci, input$strata)
  .formatResult("Fst", result, input$keep.null)
}

#' @rdname popStructStat
#' @export
#' 
statFst <- function(g, nrep = NULL, keep.null = FALSE, ...) {
  .statFst(.formatCinput(g, nrep, keep.null))
}

#-----


#' @noRd
#' 
.statFstPrime <- function(input) {
  if(is.null(input)) return(.formatResult("F'st", NULL, input$keep.null))
  if(input$ploidy == 1) return(.formatResult("F'st", NULL, input$keep.null))
  result <- statFstPrime_C(input$loci, input$strata)
  .formatResult("F'st", result, input$keep.null)
}

#' @rdname popStructStat
#' @export
#' 
statFstPrime <- function(g, nrep = NULL, keep.null = FALSE, ...) {  
  .statFstPrime(.formatCinput(g, nrep, keep.null))
}

#-----


#' @noRd
#' 
.statGst <- function(input) {
  if(is.null(input)) return(.formatResult("Gst", NULL, input$keep.null))
  if(input$ploidy == 1) return(.formatResult("Gst", NULL, input$keep.null))
  result <- statGst_C(input$loci, input$strata)
  .formatResult("Gst", result, input$keep.null)
}

#' @rdname popStructStat
#' @export
#' 
statGst <- function(g, nrep = NULL, keep.null = FALSE, ...) {  
  .statGst(.formatCinput(g, nrep, keep.null))
}

#-----


#' @noRd
#' 
.statGstPrime <- function(input) {
  if(is.null(input)) return(.formatResult("G'st", NULL, input$keep.null))
  if(input$ploidy == 1) return(.formatResult("G'st", NULL, input$keep.null))
  result <- statGstPrime_C(input$loci, input$strata, input$prime.type)
  .formatResult("G'st", result, input$keep.null)
}

#' @rdname popStructStat
#' @export
#' 
statGstPrime <- function(g, nrep = NULL, keep.null = FALSE,
                         prime.type = c("hedrick", "nei"), ...) { 
  .statGstPrime(.formatCinput(g, nrep, keep.null, ...))
}

#-----


#' @noRd
#' 
.statGstDblPrime <- function(input) {
  if(is.null(input)) return(.formatResult("G''st", NULL, input$keep.null))
  if(input$ploidy == 1) return(.formatResult("G''st", NULL, input$keep.null))
  result <- statGstDblPrime_C(input$loci, input$strata)
  .formatResult("G''st", result, input$keep.null)
}

#' @rdname popStructStat
#' @export
#' 
statGstDblPrime <- function(g, nrep = NULL, keep.null = FALSE, ...) {
  .statGstDblPrime(.formatCinput(g, nrep, keep.null))
}

#-----


#' @noRd
#' 
.statPhist <- function(input)  {
  if(is.null(input)) return(.formatResult("PHIst", NULL, input$keep.null))
  if(input$ploidy != 1) return(.formatResult("PHIst", NULL, input$keep.null))
  result <- statPhist_C(input$loci, input$strata, input$hap.dist)
  .formatResult("PHIst", result, input$keep.null)
}

#' @rdname popStructStat
#' @export
#' 
statPhist <- function(g, nrep = NULL, keep.null = FALSE,
                      model = "K80", gamma = FALSE, 
                      pairwise.deletion = TRUE, ...)  {
  .statPhist(.formatCinput(g, nrep, keep.null, hap.dist = TRUE, ...))
}

#-----

#' @rdname popStructStat
#' @export
#' 
Hstats <- function(g) {
  if(getPloidy(g) < 2) {
    result <- matrix(NA, nrow = 3, ncol = getNumLoci(g))
    rownames(result) <- c("Ho", "Hs", "Ht")
    colnames(result) <- getLociNames(g)
    return(result)
  }
  
  input <- .formatCinput(g, NULL, NULL)
  result <- Hstats_C(input$loci, input$strata)
  rownames(result) <- c("Ho", "Hs", "Ht")
  colnames(result) <- getLociNames(g)
  result
}