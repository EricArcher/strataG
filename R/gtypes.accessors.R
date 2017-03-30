#' @title \code{gtypes} Accessors
#' @description Accessors for slots in \linkS4class{gtypes} objects.
#' 
#' @param x a \linkS4class{gtypes} object.
#' @param seqName the name (or number) of a set of sequences from the 
#'   \code{@@sequences} slot to return.
#' @param as.haplotypes return sequences as haplotypes? If \code{TRUE}, contents of 
#'   \code{@@sequences} slot are returned. If \code{FALSE}, one sequence per 
#'   individual is returned.
#' @param i,j,k subsetting slots for individuals (\code{i}), loci (\code{j}),
#'   or strata (\code{k}). See Details for more information.
#' @param quiet suppress warnings about unmatched requested individuals, loci, 
#'   or strata?
#' @param drop if \code{TRUE} the return object will have unused sequences removed.
#' @param ... other arguments passed from generics (ignored).
#' @param value value being assigned by accessor.
#' 
#' @details 
#' Indexing a \code{gtypes} object with integers, characters, or logicals with 
#'   the \code{[} operator follows the same rules as normal indexing in R. The 
#'   order that individuals, loci, and strata are chosen in the order 
#'   returned by \code{indNames}, \code{locNames}, and \code{strataNames} 
#'   respectively. If unstratified samples are present, they can be selected as
#'   a group either by including \code{NA} in the character or numeric vector of the 
#'   \code{k} slot, or by providing a logical vector based on \code{is.na(strata(g))} 
#'   to the \code{i} slot.
#'
#' @return
#' \describe{
#'   \item{nInd}{number of individuals}
#'   \item{nLoc}{number of loci}
#'   \item{nStrata}{number of strata}
#'   \item{indNames}{vector of individual/sample names}
#'   \item{locNames}{vector of locus names}
#'   \item{strataNames}{vector of strata names for current scheme}
#'   \item{ploidy}{number of alleles at each locus}
#'   \item{other}{contents of \code{@@other} slot}
#'   \item{strata}{return or modify the current stratification}
#'   \item{schemes}{return or modify the current stratification schemes}
#'   \item{alleleNames}{return a list of alleles at each locus}
#'   \item{sequences}{return the \linkS4class{multidna} object in the 
#'     \code{@@sequences} slot. See \code{\link[apex]{getSequences}} to 
#'     extract individual genes or sequences from this object}
#'   \item{description}{return the object's description}
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' #--- create a diploid (microsatellite) gtypes object
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#' 
#' nStrata(msats.g)
#' strataNames(msats.g)
#' nLoc(msats.g)
#' locNames(msats.g)
#' 
#' # reassign all samples to two randomly chosen strata
#' strata(msats.g) <- sample(c("A", "B"), nInd(msats.g), rep = TRUE)
#' msats.g
#' 
#' 
#' #--- a sequence example
#' library(ape)
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' wood.g <- sequence2gtypes(x)
#' strata(wood.g) <- sample(c("A", "B"), nInd(wood.g), rep = TRUE)
#' wood.g
#' 
#' # get the multidna sequence object
#' multi.seqs <- sequences(wood.g)
#' class(multi.seqs) # "multidna"
#'
#' # get a list of DNAbin objects
#' library(apex)
#' dnabin.list <- getSequences(multi.seqs)
#' class(dnabin.list) # "list"
#' 
#' # get a DNAbin object of the first locus
#' dnabin.1 <- getSequences(multi.seqs, locNames(wood.g)[1])
#' class(dnabin.1) # "DNAbin"
#' 
#' # NOTE: The default to the 'simplify' argument in 'getSequences' is TRUE, 
#' #   so if there is only one locus, 'getSequences' will return a DNAbin object
#' #   rather than a single element list unless 'simplify = FALSE':
#' gene1 <- wood.g[, "gene1", ]
#' gene1.dnabin <- getSequences(sequences(gene1))
#' class(gene1.dnabin) # "DNAbin"
#' 
#' @name gtypes.accessors
#' @aliases accessors
#' @importFrom methods setGeneric setMethod validObject new
#' 
NULL


#' @rdname gtypes.accessors
#' @aliases nInd
#' @export
#' 
setMethod("nInd", "gtypes", function(x, ...) {
  ids <- NULL # For CRAN CHECK
  x@data[, uniqueN(ids)]
})


#' @rdname gtypes.accessors
#' @aliases nLoc
#' @export
#' 
setMethod("nLoc", "gtypes", function(x, ...) ncol(x@data) - 2)


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("nStrata", function(x, ...) standardGeneric("nStrata"))

#' @rdname gtypes.accessors
#' @aliases nStrata
#' @export
#' 
setMethod("nStrata", "gtypes", function(x, ...) x@data[, uniqueN(strata)])


#' @rdname gtypes.accessors
#' @aliases indNames
#' @export
#' 
setMethod("indNames", "gtypes", function(x, ...) {
  ids <- NULL # For CRAN CHECK
  x@data[, unique(ids)]
})


#' @rdname gtypes.accessors
#' @aliases locNames
#' @export
#' 
setMethod("locNames", "gtypes", function(x, ...) {
  setdiff(colnames(x@data), c("ids", "strata"))
})


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("strataNames", function(x, ...) standardGeneric("strataNames"))

#' @rdname gtypes.accessors
#' @aliases strataNames
#' @export
#' 
setMethod("strataNames", "gtypes", function(x, ...) {
  x@data[, sort(as.character(unique(strata)))]
})


#' @rdname gtypes.accessors
#' @aliases ploidy
#' @export
#' 
setMethod("ploidy", "gtypes", function(x, ...) x@ploidy)


#' @rdname gtypes.accessors
#' @aliases other
#' @export
#' 
setMethod("other", "gtypes", function(x, ...) x@other)


#' @rdname gtypes.accessors
#' @aliases strata
#' @export
#' 
setMethod("strata", "gtypes", function(x) {
  ids <- strata <- NULL # For CRAN CHECK
  mat <- as.matrix(x@data[, list(ids, strata)])
  mat <- mat[!duplicated(mat[, "ids"]), ]
  vec <- mat[, "strata"]
  names(vec) <- mat[, "ids"]
  vec
})


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("strata<-", function(x, value) standardGeneric("strata<-"))

#' @rdname gtypes.accessors
#' @aliases strata
#' @export
#' 
setMethod("strata<-", "gtypes", function(x, value) {
  ids <- strata <- NULL # For CRAN CHECK
  value <- if(is.null(names(value))) {
    rep(value, length.out = nrow(x@data))
  } else {
    value[x@data[, ids]]
  }
  x@data[, strata := value]
  validObject(x)
  x
})


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("schemes", function(x, ...) standardGeneric("schemes"))

#' @rdname gtypes.accessors
#' @aliases schemes
#' @export
#' 
setMethod("schemes", "gtypes", function(x, ...) x@schemes)


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("schemes<-", function(x, value) standardGeneric("schemes<-"))

#' @rdname gtypes.accessors
#' @aliases schemes
#' @export
#' 
setMethod("schemes<-", "gtypes", function(x, value) {
  x@schemes <- value
  validObject(x)
  x
})


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("alleleNames", function(x, ...) standardGeneric("alleleNames"))

#' @rdname gtypes.accessors
#' @aliases alleleNames
#' @export
#' 
setMethod("alleleNames", "gtypes", function(x) {
  sapply(x@data[, locNames(x), with = FALSE], function(x) {
    as.vector(na.omit(unique(as.character(x))))
  }, simplify = FALSE)
})


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("sequences", function(x, ...) standardGeneric("sequences"))

#' @rdname gtypes.accessors
#' @aliases sequences
#' @export
#' 
setMethod("sequences", "gtypes", function(x, seqName = NULL, as.haplotypes = TRUE, ...) {
  if(is.null(x@sequences)) return(NULL)
  dna <- getSequences(x@sequences, simplify = FALSE)
  if(!as.haplotypes) {
    dna <- lapply(locNames(x), function(l) {
      haps <- as.array(x, loci = l)
      ind.seqs <- dna[[l]][haps]
      names(ind.seqs) <- indNames(x)
      ind.seqs
    })
  }
  if(!is.null(seqName)) dna <- dna[seqName]
  dna <- as.multidna(dna)
  setLocusNames(dna) <- locNames(x)
  dna
})


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("description", function(x, ...) standardGeneric("description"))

#' @rdname gtypes.accessors
#' @aliases description
#' @export
#' 
setMethod("description", "gtypes", function(x, ...) x@description)


#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("description<-", function(x, value) standardGeneric("description<-"))

#' @rdname gtypes.accessors
#' @aliases description
#' @export
#' 
setMethod("description<-", "gtypes", function(x, value) {
  x@description <- value
  validObject(x)
  x
})


#' @rdname gtypes.accessors
#' @aliases index subset
#' @export
#' 
setMethod("[", 
          signature(x = "gtypes", i = "ANY", j = "ANY", drop = "ANY"), 
          function(x, i, j, k, ..., quiet = TRUE, drop = FALSE) {
  
  # check ids (i)
  if(missing(i)) i <- TRUE
  if(is.factor(i)) i <- as.character(i)
  ids <- indNames(x)
  i <- if(is.character(i)) {
    i <- unique(i)
    missing.ids <- setdiff(i, ids)
    if(length(missing.ids) > 0) {
      missing.ids <- paste(missing.ids, collapse = ", ")
      warning("the following ids cannot be found: ", missing.ids)
    }
    intersect(i, ids)
  } else ids[i]
  
  # check loci (j)
  if(missing(j)) j <- TRUE
  if(is.factor(j)) j <- as.character(j)
  locs <- locNames(x)
  j <- if(is.character(j)) {
    j <- unique(j)
    missing.locs <- setdiff(j, locs)
    if(length(missing.locs) > 0 & !quiet) {
      missing.locs <- paste(missing.locs, collapse = ", ")
      warning("the following loci cannot be found: ", missing.locs)
    }
    intersect(j, locs)
  } else locs[j]
  
  # check strata (k) 
  if(missing(k)) k <- TRUE
  if(is.factor(k)) k <- as.character(k)
  st <- strataNames(x)
  k <- if(is.character(k)) {
    k <- unique(k)
    missing.strata <- setdiff(k, st)
    if(length(missing.strata) > 0 & !quiet) {
      missing.strata <- paste(missing.strata, collapse = ", ")
      warning("the following strata cannot be found: ", missing.strata)
    }
    intersect(k, st)
  } else {
    st[k]
  }

  if(length(i) == 0) stop("no samples selected")
  if(length(j) == 0) stop("no loci selected")
  if(length(k) == 0) stop("no strata selected")
  
  x@data <- x@data[ids %in% i & strata %in% k, c("ids", "strata", j), with = FALSE, nomatch = 0]
  if(nrow(x@data) == 0) stop("none of the specified ids were found in the specified strata.")
  
  # check ids in selected strata
  missing.ids <- setdiff(i, x@data[, unique(ids)])
  if(length(missing.ids) > 0 & !quiet) {
    missing.ids <- paste(missing.ids, collapse = ", ")
    warning("the following ids are not in the selected strata: ", missing.ids)
  }
  
  # filter sequences
  if(!is.null(x@sequences)) {
    j.seqs <- getSequences(
      x@sequences, loci = j, simplify = FALSE, 
      exclude.gap.only = FALSE
    )
    x@sequences <- new("multidna", j.seqs)
  }
  
  # Check for samples missing data for all loci
  x <- .removeIdsMissingAllLoci(x)
  
  if(drop) x <- removeSequences(x)
  return(x)
})
