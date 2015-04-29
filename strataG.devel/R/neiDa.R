#' @title Nei's Da
#' @description Calcuate frequency-based Nei's Da for haploid or diploid data.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param na.rm logical. Delete loci which have missing data for every sample in a stratum?
#' 
#' @details Returns Nei's Da for each pair of strata.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references Nei et al 1983 Accuracy of Estimated Phylogenetic Trees from 
#'   Molecular Data. J Mol Evol 19:153-170 (eqn 7)\cr
#'   Nei, M., and S. Kumar (2000) Molecular Evolution and Phylogenetics. 
#'   Oxford University Press, Oxford. (pp. 268, eqn 13.6)
#'   
#' @export
#' 
neiDa <- function(g, na.rm = FALSE) {
  st <- strata(g)
  st.vec <- levels(st)
  if(length(st.vec) < 2) {
    stop("more than one stratum required to calculate Nei's Da") 
  }
  st.pairs <- t(combn(st.vec, 2))
  st.col <- rep(st, ploidy(g))

  Da <- apply(st.pairs, 1, function(sp) {
    loc.sum <- sapply(1:ncol(g@loci), function(i) {
      locus <- g@loci[, i]
      to.use <- st.col %in% sp & !is.na(st.col) & !is.na(locus)
      freqs <- prop.table(table(locus[to.use], droplevels(st.col[to.use])))
      sum(apply(freqs, 1, function(f) if(all(f == 0)) NA else sqrt(prod(f))))
    })
    1 - sum(loc.sum, na.rm = na.rm) / sum(!is.na(loc.sum))
  })
  result <- as.data.frame(st.pairs)
  result <- cbind(result, Da)
  colnames(result) <- c("strata.1", "strata.2", "Da")
  result
}