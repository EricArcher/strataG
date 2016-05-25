#' @title Convert Between \code{gtypes} And \code{genind} objects.
#' @description Convert a \code{gtypes} object to a \code{genind} object 
#'   and vice-versa.
#' 
#' @param x either a \linkS4class{gtypes} or \linkS4class{genind} object
#'   to convert from.
#' @param type a character string indicating the type of marker for 
#'   \linkS4class{genind} objects: 'codom' stands for 'codominant' 
#'   (e.g. microstallites, allozymes); 'PA' stands for 'presence/absence' 
#'   markers (e.g. AFLP, RAPD).
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{initialize.gtypes}, \link{df2gtypes}, 
#'   \link{sequence2gtypes}, \link{as.data.frame.gtypes}, 
#'   \link{gtypes2loci}
#' 
#' @examples
#' data(msats.g)
#' 
#' # Convert to genind
#' gi <- gtypes2genind(msats.g)
#' gi
#' 
#' # Convert to gtypes
#' gt <- genind2gtypes(gi)
#' gt
#' 
#' @name gtypes2genind
#' @export
#' 
gtypes2genind <- function(x, type = c("codom", "PA")) {
  x.mat <- as.matrix(x, one.col = TRUE, sep = "/", ids = FALSE, strata = FALSE)
  df2genind(
    X = x.mat,
    sep = "/", 
    pop = strata(x)[rownames(x.mat)],
    NA.char = NA,
    ploidy = ploidy(x),
    type = match.arg(type),
    strata = schemes(x)[rownames(x.mat), ]
  )
}


#' @rdname gtypes2genind
#' @export
#' 
genind2gtypes <- function(x) {
  gen.mat <- genind2df(x, usepop = TRUE, oneColPerAll = TRUE)
  gen.mat[gen.mat == "NA"] <- NA
  has.pop <- !is.null(x@pop)
  df2gtypes(
    x = gen.mat,
    ploidy = x@ploidy[1],
    id.col = NULL,
    strata.col = if(has.pop) 1 else NULL,
    loc.col = if(has.pop) 2 else 1,
    schemes = x@strata,
    other = other(x)
  )  
}
