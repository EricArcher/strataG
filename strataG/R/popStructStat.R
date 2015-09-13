#' @name popStructStat
#' @title Population Structure statistics
#'
#' @param g a \linkS4class{gtypes} object.
#' @param strata an optional alternative stratification vector to use in place
#'   of the stratification in \code{g}. If \code{NULL} stratification is taken 
#'   from \code{g}.
#' @param prime.type type of G'st to calculate. Can be "nei" or "hedrick".
#' @param hap.dist matrix of pairwise distances between haplotypes.
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