#' @title Convert Between \code{gtypes} And \code{loci} objects.
#' @description Convert a \code{gtypes} object to a \code{\link[pegas]{loci}} object.
#' 
#' @param x a \linkS4class{gtypes} or \code{loci} formatted object.
#' @param description a label for the \code{gtypes} object (optional).
#'  
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{initialize.gtypes}, \link{df2gtypes}, 
#'   \link{sequence2gtypes}, \link{as.data.frame.gtypes}, 
#'   \link{gtypes2genind}
#' 
#' @examples
#' data(msats.g)
#' 
#' # Convert to loci
#' lc <- gtypes2loci(msats.g)
#' lc
#' 
#' # Convert to gtypes
#' gt <- loci2gtypes(lc)
#' gt 
#' 
#' @name gtypes2loci 
#' @export
#' 
gtypes2loci <- function(x) {
  if(!is.gtypes(x)) stop("'x' must be a gtypes object")
  as.data.frame(x, one.col = TRUE, sep = "/") %>% 
    dplyr::filter(!is.na(.data$stratum)) %>% 
    dplyr::mutate_all(function(x) ifelse(is.na(x), "0/0", x)) %>% 
    dplyr::mutate_all(factor) %>% 
    tibble::column_to_rownames("id") %>% 
    as.data.frame() %>% 
    pegas::as.loci(allele.sep = "/", col.pop = 1)
}

#' @rdname gtypes2loci
#' @export
#' 
loci2gtypes <- function(x, description = NULL) {
  if(!inherits(x, "loci")) stop("'x' must be a loci object")
  mat <- as.data.frame(x) %>% 
    dplyr::mutate_all(as.character) %>% 
    dplyr::select(attr(x, "locicol")) %>% 
    alleleSplit(sep = "/")
  cbind(
    id = rownames(x),
    pop = as.character(x$population),
    mat
  ) %>% 
    df2gtypes(ploidy = 2, description = description)
}
