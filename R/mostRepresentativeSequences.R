#' @title Representative Sequences
#' @description Finds the set of sequences that represent the requested 
#'   number of clusters.
#' 
#' @param x a \code{\link[ape]{DNAbin}} object.
#' @param num.seqs number of sequences to return. If \code{NULL} (default), all
#'   sequences are returned.
#' @param model a character string specifying the evolutionary model to be used. 
#'   See \link{dist.dna} for more information.
#' @param pairwise.deletion a logical indicating whether to delete sites 
#'   with missing data. See \link{dist.dna} for more information.
#' @param as.haplotypes treat sequences as haplotypes (\code{TRUE}) or expand
#'   haplotypes to one sequence per individual (\code{FALSE}). If the latter,
#'   individual frequencies are used in cluster formation.
#' @param simplify if there is a single locus, return result in a simplified
#'   form? If \code{FALSE} a list will be returned wth one element per locus.
#' 
#' @return a vector of the sequence names.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.seqs)
#' 
#' mostRepresentativeSequences(dolph.seqs, 5)
#' 
#' mostRepresentativeSequences(dolph.seqs, 3)
#' 
#' @export
#' 
mostRepresentativeSequences <- function(
  x, num.seqs = NULL, model = "raw", pairwise.deletion = TRUE, 
  as.haplotypes = TRUE, simplify = TRUE
) { 
  
  x <- if(is.gtypes(x)) {
    getSequences(x, as.haplotypes = as.haplotypes, as.multidna = TRUE)
  } else {
    as.multidna(x)
  }
  
  result <- sapply(
    apex::getSequences(x, simplify = FALSE),
    function(dna) {
      k <- if(is.null(num.seqs)) length(dna) else min(num.seqs, length(dna))
      if(k == length(dna)) return(names(dna))
      
      # calculate distance between sequences
      seq.dist <- ape::dist.dna(
        dna, 
        model = model, 
        pairwise.deletion = pairwise.deletion, 
        as.matrix = TRUE
      )
      
      # convert distances to coordinates
      opt <- options(warn = -1)
      seq.cmd <- stats::cmdscale(seq.dist, k = length(dna) - 1)
      
      # create k clusters of sequences using k-means
      k <- min(nrow(seq.cmd), k)
      seq.cl <- stats::kmeans(seq.cmd, centers = k)$cluster
      
      ids <- tapply(names(seq.cl), seq.cl, function(ind) {
        if(length(ind) == 1) return(ind)
        if(length(ind) == 2) return(sample(ind, 1))
        # find id in this cluster that is closest to centroid
        ind.dist <- as.matrix(seq.dist[ind, ind])
        cl.cmd <- stats::cmdscale(ind.dist, k = nrow(ind.dist) - 1)
        cmd.mean <- colMeans(cl.cmd)
        dist.to.centroid <- sapply(1:nrow(cl.cmd), function(i) {
          sqrt(sum((cl.cmd[i, ] - cmd.mean) ^ 2))
        })
        ind[which.min(dist.to.centroid)]
      })
      
      options(opt)
      
      unname(ids[order(match(ids, names(dna)))])
    },
    simplify = FALSE
  )
  
  if(length(result) == 1 & simplify) result[[1]] else result
}