#' @name popStructStat
#' @title Population structure statistics
#' @description Population structure statistics
#'
#' @param g a \linkS4class{gtypes} object.
#' @param nrep number specifying number of permutation replicates to use for 
#'   permutation test.
#' @param strata.mat an optional matrix of permuted stratifications. If given, 
#'   the first column is taken to be the observed stratification and all others 
#'   are random permutations. Strata must be defined by sequential integer values 
#'   starting with 0. If \code{NULL} stratification is taken from \code{g}. This 
#'   argument is primarily used internally in \code{\link{overallTest}}.
#' @param keep.null logical. Keep the null distribution from the 
#'   permutation test?
#' @param prime.type type of G'st to calculate. Can be "nei" or "hedrick".
#' @param model,gamma,pairwise.deletion parameters passed to 
#'   \code{\link[ape]{dist.dna}}.
#' @param ... optional arguments passed to or from other functions.
#'
#' @return the named statistic estimate.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @useDynLib strataG
#' @importFrom Rcpp sourceCpp
#'
NULL


.permStrata <- function(g, nrep = NULL) {
  strata.num <- cbind(as.numeric(strata(g)) - 1)
  if(is.null(nrep)) nrep <- 0
  if(nrep < 1) return(cbind(strata.num))
  null.perm <- sapply(1:nrep, function(i) sample(strata.num))
  return(cbind(strata.num, null.perm))
}

.checkStrataMat <- function(strata.mat, g, nrep) {
  if(is.null(strata.mat)) {
    strata.mat <- .permStrata(g, nrep)
  } 
  if(!(is.matrix(strata.mat) & is.numeric(strata.mat))) {
    stop("'strata.mat' must be a numeric matrix")
  }
  if(nrow(strata.mat) != nInd(g)) {
    stop("'nrow(strata.mat)' is not equal to 'nInd(g)'")
  }
  zero.ref <- apply(strata.mat, 2, function(x) any(x == 0))
  if(!all(zero.ref)) {
    stop("'all columns in 'strata.mat' must have strata designations starting with 0")
  }
  return(strata.mat)
}

.formatResult <- function(result, stat.name, keep.null) {
  if(length(result) == 1) keep.null <- FALSE
  p.val <- if(length(result) == 1) NA else mean(result >= result[1], na.rm = TRUE) 
  return(list(
    stat.name = stat.name, 
    result = c(estimate = result[1], p.val = p.val),
    null.dist = if(keep.null) result[-1] else NULL
  ))
}