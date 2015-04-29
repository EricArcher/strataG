#' @title Nucleotide Divergence
#' @description Calculate Nei's dA between strata, and distributions of 
#'   between- and within-strata nucleotide divergence (sequence distance).
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param probs a numeric vector of probabilities of the pairwise distance 
#'   distributions with values in \code{0:1}.
#' @param ... arguments passed to \code{\link[ape]{dist.dna}} such as 
#'   \code{model} or \code{pairwise.deletion}.
#' 
#' @return a list with summaries of the \code{within} and \code{between} strata 
#'   pairwise distances including Nei's dA. 
#'   
#' @references Nei, M., and S. Kumar (2000) Molecular Evolution and 
#'   Phylogenetics. Oxford University Press, Oxford. (dA: pp. 256, eqn 12.67)
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
nucleotideDivergence <- function(g, probs = c(0, 0.025, 0.5, 0.975, 1), ...) { 
  if(is.null(g@sequences)) stop("'g' must have sequences")
  
  pair.dist.summary <- function(hap1, hap2, d) {
    pws.dist <- d[hap1, hap2]
    pws.dist <- pws.dist[lower.tri(pws.dist)]
    dist.quant <- quantile(pws.dist, probs, na.rm = TRUE)
    names(dist.quant) <- paste("pct.", probs, sep = "")
    c(mean = mean(pws.dist, na.rm = TRUE), dist.quant)
  }
  
  st <- strata(g)
  st.vec <- levels(st)
  st.pairs <- t(combn(st.vec, 2))
  hap.dist <- lapply(g@sequences@dna, dist.dna, as.matrix = TRUE, ...)
  
  result <- lapply(1:ncol(g@loci), function(i) {
    within.dist <- do.call(rbind, tapply(1:nrow(g@loci), st, function(r) {
      pair.dist.summary(g@loci[r, i], g@loci[r, i], hap.dist[[i]])
    }))
  
    between.dist <- t(apply(st.pairs, 1, function(sp) {
      btwn <- pair.dist.summary(g@loci[which(st == sp[1]), i],
                                g@loci[which(st == sp[2]), i],  
                                hap.dist[[i]]
      )
      dA <- btwn["mean"] - (sum(within.dist[sp, "mean"], na.rm = TRUE) / 2)
      c(dA = unname(dA), btwn)
    }))
    between.dist <- data.frame(st.pairs, between.dist, stringsAsFactors = FALSE)
    colnames(between.dist)[1:2] <- c("strata.1", "strata.2")  
    
    list(within = within.dist, between = between.dist) 
  })
  names(result) <- locNames(g)
  result
}