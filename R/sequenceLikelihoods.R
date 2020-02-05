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
#' @param plot generate a plot of top \code{n} most unlikely sequences.
#' @param simplify if there is a single locus, return result in a simplified
#'   form? If \code{FALSE} a list will be returned wth one element per locus.
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
#' data(dolph.haps)
#' 
#' sequenceLikelihoods(dolph.haps)
#' 
#' @export
#' 
sequenceLikelihoods <- function(x, model = "N", pairwise.deletion = FALSE, 
                                n = NULL, plot = TRUE, simplify = TRUE, ...) {
  result <- sapply(
    apex::getSequences(as.multidna(x), simplify = FALSE),
    function(dna) {
      # calculate distance between sequences
      seq.dist <- ape::dist.dna(
        dna, 
        model = model, 
        pairwise.deletion = pairwise.deletion, 
        as.matrix = TRUE
      )
      if(any(is.nan(seq.dist)) | any(is.na(seq.dist))) {
        warning("NA/NaN in pairwise distance matrix. NULL returned.")
        return(NULL)
      }
      
      dist.vec <- seq.dist[lower.tri(seq.dist)]
      dist.mean <- mean(dist.vec)
      dist.sd <- stats::sd(dist.vec)
      scale <- (dist.sd ^ 2) / dist.mean
      shape <- (dist.mean / dist.sd) ^ 2
      
      do.call(rbind, sapply(rownames(seq.dist), function(i) {
        this.dist <- seq.dist[i, ]
        this.dist <- this.dist[names(this.dist) != i]
        mean.dist <- mean(this.dist, na.rm = TRUE)
        ll <- this.dist[this.dist != 0] %>% 
          stats::dgamma(shape = shape, scale = scale) %>% 
          log() %>% 
          sum(na.rm = TRUE)
        c(mean.dist = mean.dist, neg.log.lik = -ll)
      }, simplify = FALSE)) %>% 
        tibble::as.tibble() %>% 
        dplyr::mutate(
          id = rownames(seq.dist),
          min.ll = min(.data$neg.log.lik, na.rm = TRUE),
          delta.log.lik = .data$neg.log.lik - .data$min.ll
        ) %>% 
        dplyr::arrange(dplyr::desc(.data$delta.log.lik)) %>% 
        dplyr::select(
          .data$id, .data$mean.dist, .data$neg.log.lik, .data$delta.log.lik
        ) %>% 
        as.data.frame()
    },
    simplify = FALSE
  )
  
  if(plot) {
    for(locus in names(result)) {
      df <- result[[locus]]
      if(is.null(n)) n <- nrow(df)
      n <- min(n, nrow(df))
      if(n > 0) {
        p <- df[1:n, ] %>% 
          dplyr::mutate(id = stats::reorder(.data$id, .data$delta.log.lik)) %>% 
          ggplot2::ggplot(
            ggplot2::aes_string(x = "delta.log.lik", y = "id")
          ) +
          ggplot2::geom_point() +
          ggplot2::labs(
            x = expression(paste("-", Delta, "lnL")),
            title = locus
          ) +
          ggplot2::theme(axis.title.y = ggplot2::element_blank())
        print(p)
      }
    }
  }
  
  if(length(result) == 1 & simplify) result[[1]] else result
}