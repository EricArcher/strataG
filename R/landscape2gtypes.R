#' @title Create gtypes from Rmetasim landscape
#' @description Create a gtypes object from an Rmetasim landscape object.
#' 
#' @param Rland rmetasim landscape object
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
landscape2gtypes <- function(Rland) {
  pl <- rmetasim::landscape.ploidy(Rland)
  strata <- Rland$individuals[, 1] %/% Rland$intparam$stages + 1
  gen.data <- Rland$individuals[, -(1:rmetasim::landscape.democol())]
  loc.names <- paste0("Locus", 1:length(pl))
  colnames(gen.data) <- paste(rep(loc.names, each = pl[1]), 1:pl[1], sep = ".")
  cbind(id = Rland$individuals[, 4], strata = strata, gen.data) %>% 
    df2gtypes(ploidy = pl[1])
}