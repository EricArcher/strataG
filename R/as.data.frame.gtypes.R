#' @name as.data.frame.gtypes
#' @title Convert \code{gtypes} To \code{data.frame}
#' @description Create a data.frame from a \linkS4class{gtypes} object.
#'   
#' @param x a \linkS4class{gtypes} object.
#' @param one.col logical. If \code{TRUE}, then result has one column per 
#'   locus.
#' @param sep character to use to separate alleles if \code{one.col} is 
#'   \code{TRUE}.
#' @param ids logical. include a column for individual identifiers (\code{ids})?
#' @param strata logical. include a column for current statification (\code{strata})?
#' @param sort.alleles logical. for non-haploid objects, should alleles be sorted 
#'   in genotypes or left in original order? (only takes affect if \code{one.col = TRUE})
#' @param ... additional arguments to be passed to or from methods.
#'   
#' @return A \code{data.frame} with one row per sample.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{df2gtypes} \link[strataG]{as.matrix} \link[strataG]{as.array}
#' 
#' @examples 
#' data(msats.g)
#' 
#' # with defaults (alleles in multiple columns, with ids and stratification)
#' df <- as.data.frame(msats.g)
#' str(df)
#' 
#' # one column per locus
#' onecol.df <- as.data.frame(msats.g, one.col = TRUE)
#' str(onecol.df)
#' 
#' # just the genotypes
#' genotypes.df <- as.data.frame(msats.g, one.col = TRUE, ids = FALSE, strata = FALSE)
#' str(genotypes.df)
#' 
#' @aliases as.data.frame,gtypes-method as.data.frame.gtypes as.data.frame
#' @importFrom methods setMethod
#'
#' @export
#' 
setMethod(
  "as.data.frame", "gtypes",
  function(x, one.col = FALSE, sep = "/", ids = TRUE, strata = TRUE, sort.alleles = TRUE, ...) {
    as.data.frame(as.matrix(
      x, one.col = one.col, sep = sep, ids = ids, strata = strata, sort.alleles = sort.alleles
    ), stringsAsFactors = FALSE)
  }
)
