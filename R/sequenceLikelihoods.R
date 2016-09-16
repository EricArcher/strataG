#' @title Sequence Likelihoods
#' @description Calculate likelihood of each sequence based on gamma 
#'   distribution of pairwise distances.
#' 
#' @param x a \code{\link[ape]{DNAbin}} object.
#' @param model a character string specifying the evolutionary model to be 
#'   used. Passed to \code{\link[ape]{dist.dna}}. 
#' @param pairwise.deletion a logical indicating whether to delete the 
#'   sites with missing data in a pairwise way. Passed to 
#'   \code{\link[ape]{dist.dna}}.
#' @param n number of sequences with lowest delta(log-likelihoods) to 
#'   plot. Defaults to all sequences Set to 0 to supress plotting.
#' @param ... arguments passed from other functions (ignored).
#' 
#' @details Fits a Gamma distribution to the pairwise distances of sequences 
#'   and calculates the log-likelihood for each (sum of all pairwise 
#'   log-likelihoods for that sequence). Sequences that are extremely 
#'   different from all others will have low log-likelihoods. Values returned 
#'   as delta(log-likelhoods) = difference of log-likelihoods from maximum 
#'   observed values.
#' 
#' @return vector of delta(log-Likelihoods) for each sequence, sorted from 
#'   smallest to largest, and a plot of their distributions.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' library(ape)
#' data(dolph.haps)
#' 
#' sequenceLikelihoods(as.DNAbin(dolph.haps))
#' 
#' @importFrom stats sd dgamma reorder
#' @importFrom ggplot2 ggplot aes_string geom_point xlab theme element_blank
#' @export
#' 
sequenceLikelihoods <- function(x, model = "N", pairwise.deletion = FALSE, 
                                n = NULL, ...) {
  
  if(!inherits(x, "DNAbin")) stop("'x' must be a DNAbin object")
  
  # calculate distance between sequences
  seq.dist <- dist.dna(
    x, model = model, pairwise.deletion = pairwise.deletion, as.matrix = TRUE
  )
  if(any(is.nan(seq.dist)) | any(is.na(seq.dist))) {
    warning("NA/NaN in pairwise distance matrix. NULL returned.")
    return(NULL)
  }
  
  dist.vec <- seq.dist[lower.tri(seq.dist)]
  dist.mean <- mean(dist.vec)
  dist.sd <- sd(dist.vec)
  scale <- (dist.sd ^ 2) / dist.mean
  shape <- (dist.mean / dist.sd) ^ 2
  
  mat <- t(sapply(rownames(seq.dist), function(this.seq) {
    this.dist <- seq.dist[this.seq, ]
    this.dist <- this.dist[names(this.dist) != this.seq]
    mean.dist <- mean(this.dist, na.rm = TRUE)
    this.dist <- this.dist[this.dist != 0]
    ll <- -sum(log(dgamma(this.dist, shape = shape, scale = scale)), na.rm = TRUE)
    c(mean.dist = mean.dist, negLogLik = ll)
  }))
  mat <- cbind(mat, deltaLogLik = mat[, "negLogLik"] - min(mat[, "negLogLik"], na.rm = T))
  df <- data.frame(id = rownames(mat), mat)
  df$id <- reorder(df$id, df$deltaLogLik)
  rownames(df) <- df$id
  
  if(is.null(n)) n <- nrow(df)
  n <- min(n, nrow(df))
  
  if(n > 0) {
    df.sort <- df[order(-df$deltaLogLik), ]
    p <- ggplot(df.sort[1:n, ], aes_string(x = "deltaLogLik", y = "id")) +
      geom_point() +
      xlab(expression(paste("-", Delta, "lnL"))) +
      theme(axis.title.y = element_blank())
    print(p)
  }
  
  rownames(df) <- NULL  
  df
}