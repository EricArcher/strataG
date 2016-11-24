#' @title Internal Functions
#' @description Functions intended to be used only by functions within 
#'   the strataG package.
#'   
#' @details \describe{
#'   \item{.getFileLabel}{}
#'   \item{.parseLocusNames}{}
#'   \item{.setupClusters}{}
#'   \item{.strataPairs}{}
#' }
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @name strataG-internal
NULL


#' @rdname strataG-internal
#' @param g a \linkS4class{gtypes} object.
#' @param label label for filename(s). Default is the gtypes description if present.
#' @keywords internal
#' 
.getFileLabel <- function(g, label = NULL) {
  desc <- description(g)
  label <- if(!is.null(label)) {
    label 
  } else if(!is.null(desc)) {
    desc 
  } else "strataG.gtypes"
  gsub("[[:punct:]]", ".", label)
}
  

#' @rdname strataG-internal
#' @param locus.names a vector of column names, where each locus must be
#'   named with the same roots. For example, diploid locus 'ABCD' would have
#'   two columns named something like 'ABCD.1' and 'ABCD.2', or
#'   'ABCD_A' and 'ABCD_B'.\cr
#' @param ploidy integer representing the ploidy of the data.
#' @keywords internal
#' 
.parseLocusNames <- function(locus.names, ploidy) {
  if(ploidy == 1) return(locus.names)
  loc.i <- matrix(1:length(locus.names), nrow = ploidy)
  apply(loc.i, 2, function(i) {
    this.loc <- locus.names[i]
    max.length <- max(sapply(this.loc, nchar))
    # find location of first difference
    ptr <- 1
    all.same <- length(unique(substr(this.loc, 1, ptr))) == 1
    while(all.same & ptr < max.length) {
      ptr <- ptr + 1
      all.same <- length(unique(substr(this.loc, 1, ptr))) == 1
    }
    ptr <- ptr - 1
    if(ptr == 0) {
      this.loc[1] # return first name if all names are different
    } else {
      # remove last same character if it is not alphanumeric
      if(!substr(this.loc[1], ptr, ptr) %in% c(LETTERS, letters, 0:9)) {
        ptr <- ptr - 1
      }
      substr(this.loc[1], 1, ptr)
    }
  })
}


#' @rdname strataG-internal
#' @param num.cores number of cores for multithreading. If \code{NULL}, the number 
#'   used is set to the value of \code{detectCores() - 1}.
#' @importFrom parallel detectCores makePSOCKcluster makeForkCluster
#' @keywords internal
#' 
.setupClusters <- function(num.cores = NULL) {
  # setup clusters
  if(is.null(num.cores)) num.cores <- detectCores() - 1
  if(is.na(num.cores)) num.cores <- 1
  num.cores <- max(1, num.cores)
  num.cores <- min(num.cores, detectCores() - 1)
  if(num.cores > 1) {
    is.windows <- .Platform$OS.type == "windows"
    cl.func <- ifelse(is.windows, makePSOCKcluster, makeForkCluster)
    cl.func(num.cores)
  } else NULL
}


#' @rdname strataG-internal
#' @param g a \linkS4class{gtypes} object.
#' @importFrom utils combn
#' @keywords internal
#' 
.strataPairs <- function(g) {
  st <- strataNames(g)
  if(length(st) < 2) return(NULL)
  strata.pairs <- t(combn(st, 2))
  colnames(strata.pairs) <- c("strata.1", "strata.2")
  as.data.frame(strata.pairs, stringsAsFactors = FALSE)
}


#' @rdname strataG-internal
#' @param g a \linkS4class{gtypes} object.
#' @keywords internal
#' 
.removeIdsMissingAllLoci <- function(g) {
  ids <- NULL # For CRAN CHECK
  is.missing.all <- g@data[, list(missing = all(is.na(.SD))), .SDcols = !c("ids", "strata"), by = "ids"]
  to.remove <- is.missing.all[(missing)]$ids
  if(length(to.remove) > 0) {
    warning(
      "The following samples are missing data for all loci and have been removed: ", 
      paste(to.remove, collapse = ", ")
    )
    g@data <- g@data[!ids %in% to.remove]
  }
  g@data <- droplevels(g@data)
  g
}


#' @rdname strataG-internal
#' @param fun a function that takes one locus column at a time.
#' @param g a \linkS4class{gtypes} object.
#' @keywords internal
#' 
.applyPerLocus <- function(fun, g, ...) {
  g@data[, apply(.SD, 2, fun, ...), .SDcols = !c("ids", "strata")]
}


#' @rdname strataG-internal
#' @param g a \linkS4class{gtypes} object.
#' @keywords internal
#' 
.numericLoci <- function(g, min.val = 0) {
  ids <- NULL # For CRAN CHECK
  .convToNum <- function(x) min.val + (as.numeric(droplevels(x)) - 1)
  list(
    ids = g@data[, unique(ids)],
    loci = as.matrix(g@data[, sapply(.SD, .convToNum), .SDcols = !c("ids", "strata")])
  )
}