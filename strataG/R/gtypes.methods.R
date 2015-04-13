setMethod("nInd", "gtypes", function(x) nrow(x@loci) / x@ploidy)
setMethod("nLoc", "gtypes", function(x) ncol(x@loci))
setMethod("indNames", "gtypes", function(x) {
  ids <- rownames(x@loci)[1:(nrow(x@loci) / x@ploidy)]
  substr(ids, 1, nchar(ids) - 2)
})
setMethod("locNames", "gtypes", function(x) colnames(x@loci))
setMethod("ploidy", "gtypes", function(x) x@ploidy)
setMethod("other", "gtypes", function(x) x@other)

setGeneric("strata", function(g, ...) standardGeneric("strata"))
setMethod("strata", "gtypes", function(g) g@strata)

setGeneric("schemes", function(g, ...) standardGeneric("schemes"))
setMethod("schemes", "gtypes", function(g) g@schemes)

setGeneric("schemes<-", function(g, value) standardGeneric("schemes<-"))
setMethod("schemes<-", "gtypes", function(g, value) {
  g@schemes <- value
  validObject(g)
  g
})
