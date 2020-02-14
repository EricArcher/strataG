#' @name landscape2gtypes
#' 
#' @title Convert Rmetasim landscape
#' @description `landscape2gtypes` creates a gtypes object from an 
#'   Rmetasim landscape object. `landscape2df` creates a data.frame.
#' 
#' @param Rland rmetasim landscape object
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
landscape2gtypes <- function(Rland) {
  df2gtypes(landscape2df(Rland), ploidy = rmetasim::landscape.ploidy(Rland)[1])
}

#' @rdname landscape2gtypes
#' @export
#' 
landscape2df <- function(Rland) {
  pl <- rmetasim::landscape.ploidy(Rland)
  strata <- Rland$individuals[, 1] %/% Rland$intparam$stages + 1
  gen.data <- Rland$individuals[, -(1:rmetasim::landscape.democol()), drop = FALSE]
  loc.names <- paste0("Locus", 1:length(pl))
  colnames(gen.data) <- paste(rep(loc.names, each = pl[1]), 1:pl[1], sep = ".")
  cbind(id = Rland$individuals[, 4], strata = strata, gen.data) %>% 
    as.data.frame(stringsAsFactors = FALSE)
}