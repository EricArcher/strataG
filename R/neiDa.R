#' @title Nei's Da
#' @description Calcuate frequency-based Nei's Da for haploid or diploid data.
#' 
#' @param g a \linkS4class{gtypes} object.
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
#' @examples
#' data(msats.g)
#' 
#' neiDa(msats.g)
#'  
#' @export
#' 
neiDa <- function(g) {
  if(getNumStrata(g) < 2) stop("cannot compute Nei's Da with < 2 strata")
  g <- .checkHapsLabelled(g)
  
  freqs <- alleleFreqs(g, by.strata = TRUE)
  .DaFunc <- function(st1, st2, freqs) {    
    loc.sum <- purrr::map_dbl(freqs, function(loc.freqs) {
      freqs.st <- prop.table(loc.freqs[, c(st1, st2), drop = FALSE])
      apply(freqs.st, 1, function(f) if(all(f == 0)) NA else sqrt(prod(f))) %>% 
        sum(na.rm = TRUE)
    }) 
    1 - sum(loc.sum, na.rm = TRUE) / sum(!is.na(loc.sum))
  }
  
  .strataPairs(g) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(Da = .DaFunc(.data$strata.1, .data$strata.2, freqs)) %>% 
    as.data.frame()
}