#' @title Population structure statistics
#' @description Population structure statistics
#'
#' @param g a \linkS4class{gtypes} object.
#' @param nrep number specifying number of permutation replicates to use for 
#'   permutation test.
#' @param strata.mat an optional matrix of permuted stratifications. See Notes 
#'  for more details. Ignored if \code{nrep} is not \code{NULL}. 
#' @param keep.null logical. Keep the null distribution from the 
#'   permutation test?
#' @param prime.type type of G'st to calculate. Can be "nei" or "hedrick".
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
#'   \item{Hstats}{Nei, M. and R.K. Chesser. 1983. Estimation of fixation 
#'     indices and gene diversities. Ann. Hum. Genet. 47:253-259.}
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @name popStructStat
#' @useDynLib strataG, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
NULL

#----- Supporting Functions -----

# create matrix of permuted stratifications
.permStrata <- function(g, nrep = NULL) {
  ids <- strata <- NULL # For CRAN CHECK
  st <- g@data[, list(st = unique(strata)), by = ids]$st
  if(any(is.na(st))) stop("cannot run with unstratified samples")
  strata.num <- cbind(as.numeric(factor(st)) - 1)
  rownames(strata.num) <- indNames(g)
  if(is.null(nrep)) nrep <- 0
  if(nrep < 1) return(cbind(strata.num))
  st.mat <- cbind(strata.num, sapply(1:nrep, function(i) sample(strata.num)))
  rownames(st.mat) <- rownames(strata.num)
  st.mat
}

# check that matrix of permuted stratifications is of correct format
.checkStrataMat <- function(strata.mat, g, nrep = NULL) {
  if(is.null(strata.mat)) {
    strata.mat <- .permStrata(g, nrep)
  } else {
    if(!(is.matrix(strata.mat) & is.numeric(strata.mat))) {
      stop("'strata.mat' must be a numeric matrix")
    }
    if(nrow(strata.mat) != nInd(g)) {
      stop("'nrow(strata.mat)' is not equal to 'nInd(g)'")
    }
    for(i in 1:ncol(strata.mat)) {
      if(!any(strata.mat[, i] == 0)) {
        stop("column ", i, " in 'strata.mat' does not have strata designations starting with 0")
      }
      if(any(is.na(strata.mat[, i]))) {
        stop("column ", i, " has unstratified samples")
      }
    }
  }
  return(strata.mat)
}

# create a list of formatted numeric loci and strata matrices for input to C
#   population structure functions
.formatCinput <- function(g, strata.mat, nrep) {
  ids <- NULL # For CRAN CHECK
  setkey(g@data, ids)
  loci <- .numericLoci(g)
  strata <- .checkStrataMat(strata.mat, g, nrep)
  strata <- strata[loci$ids, , drop = FALSE]
  list(loci = loci$loci, strata = strata)
}

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

#' @rdname popStructStat
#' @export
#' 
Hstats <- function(g) {
  if(ploidy(g) < 2) {
    result <- matrix(NA, nrow = 3, ncol = nLoc(g))
    rownames(result) <- c("Ho", "Hs", "Ht")
    colnames(result) <- locNames(g)
    return(result)
  }
  
  input <- .formatCinput(g, NULL, NULL)
  result <- Hstats_C(input$loci, input$strata)
  rownames(result) <- c("Ho", "Hs", "Ht")
  colnames(result) <- locNames(g)
  result
}


#' @rdname popStructStat
#' @export
#' 
statChi2 <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE, ...) {
  input <- .formatCinput(g, strata.mat, nrep)
  result <- statChi2_C(input$loci, input$strata)
  return(.formatResult("Chi2", result, keep.null))
}


#' @rdname popStructStat
#' @export
#' 
statFis <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE, ...) {
  if(ploidy(g) == 1 | nStrata(g) == 1) {
    return(.formatResult("Fis", NULL, keep.null))
  }
  
  input <- .formatCinput(g, strata.mat, nrep)
  result <- statFstPrime_C(input$loci, input$strata)
  return(.formatResult("Fis", result, keep.null))
}


#' @rdname popStructStat
#' @export
#' 
statFst <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE, ...) {
  # for haploid data, run Phist with no sequences
  if(ploidy(g) == 1) {
    g@sequences <- NULL
    result <- statPhist(
      g, nrep = nrep, strata.mat = strata.mat, keep.null = keep.null
    )
    result$stat.name <- "Fst"
    return(result)
  }
  
  input <- .formatCinput(g, strata.mat, nrep)
  result <- statFst_C(input$loci, input$strata)
  return(.formatResult("Fst", result, keep.null))
}


#' @rdname popStructStat
#' @export
#' 
statFstPrime <- function(g, nrep = NULL, strata.mat = NULL, 
                         keep.null = FALSE, ...) {  
  if(ploidy(g) == 1 | nStrata(g) == 1) {
    return(.formatResult("F'st", NULL, keep.null))
  }
  
  input <- .formatCinput(g, strata.mat, nrep)
  result <- statFstPrime_C(input$loci, input$strata)
  return(.formatResult("F'st", result, keep.null))
}


#' @rdname popStructStat
#' @export
#' 
statGst <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE, ...) {
  if(ploidy(g) == 1 | nStrata(g) == 1) {
    return(.formatResult("Gst", NULL, keep.null))
  }
  
  input <- .formatCinput(g, strata.mat, nrep)
  result <- statGst_C(input$loci, input$strata)
  return(.formatResult("Gst", result, keep.null))
}


#' @rdname popStructStat
#' @export
#' 
statGstPrime <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE,
                         prime.type = c("nei", "hedrick"), ...) { 
  if(ploidy(g) == 1 | nStrata(g) == 1) {
    return(.formatResult("G'st", NULL, keep.null))
  }
  
  input <- .formatCinput(g, strata.mat, nrep)
  prime.type <- switch(match.arg(prime.type), nei = 0, hedrick = 1)
  result <- statGstPrime_C(input$loci, input$strata, prime.type)
  return(.formatResult("G'st", result, keep.null))
}


#' @rdname popStructStat
#' @export
#' 
statGstDblPrime <- function(g, nrep = NULL, strata.mat = NULL, 
                            keep.null = FALSE, ...) {
  if(ploidy(g) == 1 | nStrata(g) == 1) {
    return(.formatResult("G''st", NULL, keep.null))
  }
  
  input <- .formatCinput(g, strata.mat, nrep)
  result <- statGstDblPrime_C(input$loci, input$strata)
  return(.formatResult("G''st", result, keep.null))
}


#' @rdname popStructStat
#' @export
#' 
statJostD <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE, ...) {
  if(ploidy(g) == 1 | nStrata(g) == 1) {
    return(.formatResult("D", NULL, keep.null))
  }
  
  input <- .formatCinput(g, strata.mat, nrep)
  result <- statJostD_C(input$loci, input$strata)
  return(.formatResult("D", result, keep.null))
}


#' @rdname popStructStat
#' @export
#' 
statPhist <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE,
                      model = "K80", gamma = FALSE, pairwise.deletion = TRUE, 
                      ...)  {
  if(ploidy(g) != 1 | nStrata(g) == 1) {
    return(.formatResult("PHIst", NULL, keep.null))
  }
  
  hap.dist <- if(is.null(sequences(g))) {
    # format distances for Fst (all 1s, and 0s on the diagonal)
    lapply(locNames(g), function(gene) {
      haps <- alleleNames(g)[[gene]]
      hd <- matrix(1, nrow = length(haps), ncol = length(haps), 
                   dimnames = list(haps, haps))
      diag(hd) <- 0
      hd
    })
  } else {
    lapply(locNames(g), function(gene) {
      haps <- alleleNames(g)[[gene]]
      dna <- as.list(getSequences(sequences(g, gene), simplify = TRUE))[haps]
      hd <- dist.dna(
        dna, model = model, gamma = gamma, 
        pairwise.deletion = pairwise.deletion, as.matrix = TRUE
      )
      hd[haps, haps]
    })
  }
  names(hap.dist) <- locNames(g)
  
  input <- .formatCinput(g, strata.mat, nrep)
  result <- statPhist_C(input$loci, input$strata, hap.dist)
  return(.formatResult("PHIst", result, keep.null))
}