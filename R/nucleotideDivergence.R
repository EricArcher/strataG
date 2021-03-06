#' @title Nucleotide Divergence
#' @description Calculate distributions of between- and within-strata 
#'   nucleotide divergence (sequence distance), which includes 
#'   Nei's \eqn{\pi} (usually referred to as "nucleotide diversity") and 
#'   Nei's dA between strata. 
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param probs a numeric vector of probabilities of the pairwise distance 
#'   distributions with values in \code{0:1}.
#' @param model evolutionary model to be used. see \code{\link[ape]{dist.dna}} 
#'   for options.
#' @param ... other arguments passed to \code{\link[ape]{dist.dna}}.
#' 
#' @return a list with summaries of the \code{$within} and \code{$between} strata 
#'   pairwise distances including Nei's dA (in \code{$between}). 
#'   Nei's \eqn{\pi} is the mean between-strata divergence.
#'   
#' @references Nei, M., and S. Kumar (2000) Molecular Evolution and 
#'   Phylogenetics. Oxford University Press, Oxford. (dA: pp. 256, eqn 12.67)
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dloop.g)
#' 
#' nd <- nucleotideDivergence(dloop.g)
#' nd$within
#' nd$between
#' 
#' @aliases dA
#' @export
#' 
nucleotideDivergence <- function(g, probs = c(0, 0.025, 0.5, 0.975, 1), 
                                 model = "raw", ...) { 
  if(getPloidy(g) > 1) stop("'g' must be haploid")
  if(is.null(g@sequences)) stop("'g' must have sequences")
  
  .pair.dist.smry <- function(haps, d, probs) {
    pws.dist <- apply(haps, 2, function(h) {
      if(any(is.na(h))) return(NA)
      d[h[1], h[2]]
    })
    dist.quant <- stats::quantile(pws.dist, probs, na.rm = TRUE)
    names(dist.quant) <- paste0("q.", probs)
    c(mean = mean(pws.dist, na.rm = TRUE), dist.quant)
  }
  
  .expand.smry.cols <- function(df) {
    dplyr::bind_cols(
      df, 
      as.data.frame(do.call(rbind, df$smry))
    ) %>% 
      dplyr::select(-.data$smry) %>% 
      as.data.frame()
  }
  
  g <- g[, , getStrataNames(g)]
  hap.dist <- sapply(
    getSequences(g, as.haplotypes = TRUE, simplify = FALSE), 
    ape::dist.dna, 
    model = model, 
    as.matrix = TRUE, 
    ...,
    simplify = FALSE
  )
  bad.dists <- sapply(hap.dist, function(x) any(is.nan(x) | is.infinite(x)))
  if(any(bad.dists)) {
    warning(
      "NaN or InF returned for some pairwise distances. ",
      "See Note in ?ape::dist.dna for an explanation."
    )
  }
  
  within <- g@data %>% 
    dplyr::group_by(.data$locus, .data$stratum) %>% 
    dplyr::do(smry = {
      if(nrow(.data) == 1) {
        as.numeric(c(
          mean = NA,
          stats::setNames(rep(NA, length(probs)), paste0("q.", probs))
        ))
      } else {
        haps <- utils::combn(.data$allele, 2)
        loc <- unique(.data$locus)
        .pair.dist.smry(haps, hap.dist[[loc]], probs)
      }
    }) 
  within <- .expand.smry.cols(within)
  
  st.pairs <- .strataPairs(g)
  between <- if(is.null(st.pairs)) NA else {
    result <- do.call(rbind, lapply(getLociNames(g), function(loc) {
      cbind(locus = loc, st.pairs)
    })) %>% 
      dplyr::group_by(.data$locus, .data$strata.1, .data$strata.2) %>% 
      dplyr::do(smry = {
        st.1 <- .data$strata.1
        st.2 <- .data$strata.2
        loc <- unique(.data$locus)
        h1 <- g@data %>% 
          dplyr::filter(.data$stratum == st.1) %>% 
          dplyr::pull(.data$allele)
        h2 <- g@data %>% 
          dplyr::filter(.data$stratum == st.2) %>% 
          dplyr::pull(.data$allele)
        haps <- t(expand.grid(h1, h2))
        smry <- .pair.dist.smry(haps, hap.dist[[loc]], probs)
        wthn.sum <- within %>%
          dplyr::filter(.data$locus == loc & .data$stratum %in% c(st.1, st.2)) %>%
          dplyr::summarize(sum = sum(.data$mean, na.rm = TRUE)) %>%
          dplyr::pull("sum")
        dA <- unname(smry["mean"] - (wthn.sum / 2))
        c(dA = dA, smry)
      })
    .expand.smry.cols(result)
  }
  
  list(within = within, between = between)
}