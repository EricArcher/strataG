#' @name gtypes.accessors
#' @title gtypes Accessors
#' @description Accessors for slots in \linkS4class{gtypes} objects.
#' 
#' @param x a \linkS4class{gtypes} object.
#' @param seqName the name (or number) of a set of sequences from the 
#'   \code{@@sequences} slot to return.
#' @param ... other arguments passed from generics (ignored).
#' @param value value being assigned by accessor.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
setClass("gtypes")

#' @aliases nInd,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("nInd", "gtypes", function(x, ...) nrow(x@loci) / x@ploidy)

#' @aliases nLoc,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("nLoc", "gtypes", function(x, ...) ncol(x@loci))

#' @rdname gtypes.accessors
#' @export
setGeneric("nStrata", function(x, ...) standardGeneric("nStrata"))
#' @aliases nStrata,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("nStrata", "gtypes", function(x, ...) nlevels(x@strata))

#' @aliases indNames,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("indNames", "gtypes", function(x, ...) {
  ids <- rownames(x@loci)[1:(nrow(x@loci) / x@ploidy)]
  unique(substr(ids, 1, nchar(ids) - 2))
})

#' @aliases locNames,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("locNames", "gtypes", function(x, ...) colnames(x@loci))

#' @rdname gtypes.accessors
#' @export
setGeneric("strataNames", function(x, ...) standardGeneric("strataNames"))
#' @aliases strataNames,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("strataNames", "gtypes", function(x, ...) levels(x@strata))

#' @aliases ploidy,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("ploidy", "gtypes", function(x, ...) x@ploidy)

#' @aliases other,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("other", "gtypes", function(x, ...) x@other)

#' @rdname gtypes.accessors
#' @export
setGeneric("strata", function(x, ...) standardGeneric("strata"))
#' @aliases strata,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("strata", "gtypes", function(x, ...) x@strata)

#' @rdname gtypes.accessors
#' @export
setGeneric("strata<-", function(x, value) standardGeneric("strata<-"))
#' @aliases strata<-,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("strata<-", "gtypes", function(x, value) {
  strata <- factor(rep(value, length.out = nInd(x)))
  names(strata) <- indNames(x)
  x@strata <- droplevels(strata)
  validObject(x)
  x
})

#' @rdname gtypes.accessors
#' @export
setGeneric("schemes", function(x, ...) standardGeneric("schemes"))
#' @aliases schemes,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("schemes", "gtypes", function(x, ...) x@schemes)

#' @rdname gtypes.accessors
#' @export
setGeneric("schemes<-", function(x, value) standardGeneric("schemes<-"))
#' @aliases schemes<-,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("schemes<-", "gtypes", function(x, value) {
  x@schemes <- value
  validObject(x)
  x
})

#' @rdname gtypes.accessors
#' @export
setGeneric("seqNames", function(x, ...) standardGeneric("seqNames"))
#' @aliases seqNames,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("seqNames", "gtypes", function(x, ...) names(x@sequences@dna))

#' @rdname gtypes.accessors
#' @export
setGeneric("sequences", function(x, ...) standardGeneric("sequences"))
#' @aliases sequences,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("sequences", "gtypes", function(x, seqName = NULL, ...) {
  if(is.null(seqName)) {
    x@sequences
  } else {
    x@sequences@dna[[seqName]]
  }
})

#' @rdname gtypes.accessors
#' @export
setGeneric("description", function(x, ...) standardGeneric("description"))
#' @aliases description,gtypes
#' @rdname gtypes.accessors
#' @export
setMethod("description", "gtypes", function(x, ...) x@description)
