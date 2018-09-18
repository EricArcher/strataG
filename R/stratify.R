#' @title Stratify gtypes
#' @description Choose a new stratification scheme from the \code{schemes}
#'   slot in a \linkS4class{gtypes} object.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param scheme either the column name of a stratification scheme stored 
#'   in the data.frame of the \code{schemes} slot of \code{g}, or a vector or 
#'   factor identifying which stratum each sample belongs to.
#' @param drop remove samples not assigned to a stratum? (those assigned \code{NA} 
#'   in stratification scheme)
#'
#' @note If \code{scheme} is a vector or factor and has names, then the 
#'   they will be used to match with \code{\link{indNames}} of \code{g}. 
#'   Otherwise \code{scheme} should be the same length as the number of 
#'   samples in \code{g} or values in \code{scheme} will be recycled as 
#'   necessary.
#'   
#' @return A new \linkS4class{gtypes} object with an updated \code{strata} slot.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{schemes}
#'
#' @examples
#' data(msats.g)
#' msats.g
#' 
#' broad.msats <- stratify(msats.g, "broad")
#' broad.msats
#' 
#' @export
#' 
stratify <- function(g, scheme = NULL, drop = TRUE) {
  ids <- getIndNames(g)
  
  scheme <- if(is.null(scheme)) {
    tibble(id = ids, .new = rep("Default", length(ids)))
  } else if(!(is.vector(scheme) | is.factor(scheme))) {
    stop("'scheme' must be a vector or a factor")
  } else if(length(scheme) != 1) {
    stop("'scheme' must be one element long")
  } else if(!scheme %in% colnames(schemes(g))) {
    stop(paste("scheme '", scheme, "' cannot be found", sep = ""))
  } else {
    schemes(g) %>% 
      rename('.new' = scheme) %>% 
      select(id, .new)
  }
      
  g@data <- g@data %>% 
    left_join(scheme, by = "id") %>% 
    select(id, .new, locus, allele) %>% 
    rename(stratum = .new) %>% 
    as.data.table()

  if(drop) {
    g@data <- g@data %>% 
      filter(!is.na(stratum)) %>% 
      as.data.table()
    g <- removeSequences(g)
  } 
  g
}