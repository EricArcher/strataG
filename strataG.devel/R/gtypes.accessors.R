#' @name gtypes.accessors
#' @title \code{gtypes} Accessors
#' @description Accessors for slots in \linkS4class{gtypes} objects.
#' 
#' @param x a \linkS4class{gtypes} object.
#' @param seqName the name (or number) of a set of sequences from the 
#'   \code{@@sequences} slot to return.
#' @param ids vector of individual ids.
#' @param loci vecor of loci.
#' @param ... other arguments passed from generics (ignored).
#' @param value value being assigned by accessor.
#' 
#' @return
#' \tabular{ll}{
#'   \code{nInd} \tab number of individuals/samples.\cr
#'   \code{nLoc} \tab number of loci.\cr
#'   \code{nStrata} \tab number of strata.\cr
#'   \code{indNames} \tab vector of individual/sample names.\cr
#'   \code{locNames} \tab vector of locus names.\cr
#'   \code{strataNames} \tab vector of strata names for current scheme.\cr
#'   \code{ploidy} \tab number of alleles at each locus.\cr
#'   \code{other} \tab contents of \code{@@other} slot.\cr
#'   \code{strata} \tab return or modify the current stratification.\cr
#'   \code{schemes} \tab return or modify the current stratification schemes.\cr
#'   \code{loci} \tab return a data.frame of the alleles for the specified ids 
#'     and loci.\cr
#'   \code{seqNames} \tab return the names of each set of sequences.\cr
#'   \code{sequences} \tab return the \linkS4class{multidna} object in 
#'     the \code{@@sequences} slot.\cr
#'   \code{description} \tab return the object's description.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @aliases accessors
#' 
setClass("gtypes")

#' @rdname gtypes.accessors
#' @aliases nInd,gtypes
#' @export
setMethod("nInd", "gtypes", function(x, ...) nrow(x@loci) / x@ploidy)

#' @rdname gtypes.accessors
#' @aliases nLoc,gtypes
#' @export
setMethod("nLoc", "gtypes", function(x, ...) ncol(x@loci))

#' @rdname gtypes.accessors
#' @export
setGeneric("nStrata", function(x, ...) standardGeneric("nStrata"))
#' @rdname gtypes.accessors
#' @aliases nStrata,gtypes
#' @export
setMethod("nStrata", "gtypes", function(x, ...) nlevels(x@strata))

#' @rdname gtypes.accessors
#' @aliases indNames,gtypes
#' @export
setMethod("indNames", "gtypes", function(x, ...) {
  ids <- rownames(x@loci)[1:(nrow(x@loci) / x@ploidy)]
  unique(substr(ids, 1, nchar(ids) - 2))
})

#' @rdname gtypes.accessors
#' @aliases locNames,gtypes
#' @export
setMethod("locNames", "gtypes", function(x, ...) colnames(x@loci))

#' @rdname gtypes.accessors
#' @export
setGeneric("strataNames", function(x, ...) standardGeneric("strataNames"))
#' @rdname gtypes.accessors
#' @aliases strataNames,gtypes
#' @export
setMethod("strataNames", "gtypes", function(x, ...) levels(x@strata))

#' @rdname gtypes.accessors
#' @aliases ploidy,gtypes
#' @export
setMethod("ploidy", "gtypes", function(x, ...) x@ploidy)

#' @rdname gtypes.accessors
#' @aliases other,gtypes
#' @export
setMethod("other", "gtypes", function(x, ...) x@other)

#' @rdname gtypes.accessors
#' @export
setGeneric("strata", function(x, ...) standardGeneric("strata"))
#' @rdname gtypes.accessors
#' @aliases strata,gtypes
#' @export
setMethod("strata", "gtypes", function(x, ...) x@strata)

#' @rdname gtypes.accessors
#' @export
setGeneric("strata<-", function(x, value) standardGeneric("strata<-"))
#' @rdname gtypes.accessors
#' @aliases strata<-,gtypes
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
#' @rdname gtypes.accessors
#' @aliases schemes,gtypes
#' @export
setMethod("schemes", "gtypes", function(x, ...) x@schemes)

#' @rdname gtypes.accessors
#' @export
setGeneric("schemes<-", function(x, value) standardGeneric("schemes<-"))
#' @rdname gtypes.accessors
#' @aliases schemes<-,gtypes
#' @export
setMethod("schemes<-", "gtypes", function(x, value) {
  x@schemes <- value
  validObject(x)
  x
})

#' @rdname gtypes.accessors
#' @export
setGeneric("loci", function(x, ...) standardGeneric("loci"))
#' @rdname gtypes.accessors
#' @aliases loci,gtypes
#' @export
setMethod("loci", "gtypes", function(x, ids = NULL, loci = NULL) {
  if(is.null(ids)) ids <- indNames(x)
  if(is.null(loci)) loci <- locNames(x)
  if(!all(ids %in% indNames(x))) stop("some 'ids' not found in 'x'")
  if(!all(loci %in% locNames(x))) stop("some 'loci' not found in 'x'")
  x@loci[idRows(ids, rownames(x@loci)), loci, drop = FALSE]
})

#' @rdname gtypes.accessors
#' @export
setGeneric("seqNames", function(x, ...) standardGeneric("seqNames"))
#' @rdname gtypes.accessors
#' @aliases seqNames,gtypes
#' @export
setMethod("seqNames", "gtypes", function(x, ...) names(x@sequences@dna))

#' @rdname gtypes.accessors
#' @export
setGeneric("sequences", function(x, ...) standardGeneric("sequences"))
#' @rdname gtypes.accessors
#' @aliases sequences,gtypes
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
#' @rdname gtypes.accessors
#' @aliases description,gtypes
#' @export
setMethod("description", "gtypes", function(x, ...) x@description)
