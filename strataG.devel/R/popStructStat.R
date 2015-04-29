#' @name popStructStat
#' @title Population Structure statistics
#'
#' @param g a \linkS4class{gtypes} object.
#' @param strata an optional alternative stratification vector to use in place
#'   of the stratification in \code{g}.
#' @param prime.type type of G'st to calculate. Can be "nei" or "hedrick".
#' @param hap.dist matrix of pairwise distances between haplotypes.
#' @param pairwise.deletion logical. Do pairwise deletion of sites with 
#'   missing data. See \code{\link[ape]{dist.dna}}.
#' @param ... optional arguments passed to or from other functions.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
NULL