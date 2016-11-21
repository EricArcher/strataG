#' @title Convert \code{gtypes} To \code{matrix}
#' @description Create a matrix from a \linkS4class{gtypes} object.
#'   
#' @param x a \linkS4class{gtypes} object.
#' @param one.col logical. If \code{TRUE}, then result has one column per 
#'   locus.
#' @param sep character to use to separate alleles if \code{one.col} is 
#'   \code{TRUE}.
#' @param ids logical. include a column for individual identifiers (\code{ids})?
#' @param strata logical. include a column for current statification (\code{strata})?
#' @param sort.alleles logical. for non-haploid objects, should alleles be sorted 
#'   in genotypes or left as in original order? (only takes affect if \code{one.col = TRUE})
#' @param ... additional arguments to be passed to or from methods.
#'   
#' @return A \code{matrix} with one row per sample.
#' 
#' @author Eric Archer \email{eric.archer@@noa.gov}
#' 
#' @seealso \link{df2gtypes} \link[strataG]{as.data.frame} \link[strataG]{as.array}
#' 
#' @examples 
#' data(msats.g)
#' 
#' # with defaults (alleles in multiple columns, with ids and stratification)
#' mat <- as.matrix(msats.g)
#' head(mat)
#' 
#' # one column per locus
#' onecol.mat <- as.matrix(msats.g, one.col = TRUE)
#' head(onecol.mat)
#' 
#' # just the genotypes
#' genotypes.mat <- as.matrix(msats.g, one.col = TRUE, ids = FALSE, strata = FALSE)
#' head(genotypes.mat)
#' 
#' @aliases as.matrix,gtypes-method as.matrix.gtypes as.matrix
#' @importFrom methods setMethod
#' 
#' @export
#' 
setMethod(
  "as.matrix", "gtypes",
  function(x, one.col = FALSE, sep = "/", ids = TRUE, strata = TRUE, sort.alleles = TRUE, ...) {
    setkey(x@data, ids)
    # loop through each locus
    gen.mat <- do.call(cbind, lapply(locNames(x), function(locus) {
      # collapse alleles to create a one column matrix if one.col == TRUE
      this.loc <- if(one.col | ploidy(x) == 1) {
        arr <- as.array(x, loci = locus, drop = FALSE)
        cbind(apply(arr, 1, function(alleles) {
          if(any(is.na(alleles))) NA else {
            if(sort.alleles) alleles <- sort(alleles)
            paste(alleles, collapse = sep)
          }
        }))
      } else {
        as.array(x, loci = locus, drop = TRUE)
      }
      # assign column names
      colnames(this.loc) <- if(dim(this.loc)[2] == 1) {
        locus
      } else {
        paste(locus, 1:ncol(this.loc), sep = ".")
      } 
      this.loc
    }))
    rownames(gen.mat) <- indNames(x)
    if(strata) gen.mat <- cbind(strata = as.character(strata(x)), gen.mat)
    if(ids) gen.mat <- cbind(ids = indNames(x), gen.mat)
    gen.mat
  }
)