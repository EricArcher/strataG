#' @title \code{gtypes} Accessors
#' @description Accessors for slots in \linkS4class{gtypes} objects.
#' 
#' @param x a \linkS4class{gtypes} object.
#' @param by.strata logical - return results by strata?
#' @param seqName the name (or number) of a set of sequences from the 
#'   \code{@@sequences} slot to return.
#' @param as.haplotypes return sequences as haplotypes? If \code{TRUE}, contents of 
#'   \code{@@sequences} slot are returned. If \code{FALSE}, one sequence per 
#'   individual is returned.
#' @param as.multidna return sequences as a \linkS4class{multidna} object? If 
#'   \code{FALSE}, sequences are returned as a list.
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
#' @aliases getNumInd
#' @export
#' 
setMethod("getNumInd", "gtypes", function(x, by.strata = FALSE, ...) {
  if(by.strata) {
    x@data %>% 
      dplyr::group_by(stratum) %>% 
      dplyr::summarize(num.ind = n_distinct(id)) %>% 
      dplyr::ungroup() %>% 
      as.data.frame()
  } else {
    length(getIndNames(x))
  }
})

#' @rdname gtypes.accessors
#' @aliases getNumLoci
#' @export
#' 
setMethod("getNumLoci", "gtypes", function(x, ...) {
  x@data$locus %>% 
    dplyr::n_distinct(na.rm = TRUE)
})

#' @rdname gtypes.accessors
#' @aliases getNumStrata
#' @export
#' 
setGeneric("getNumStrata", function(x, ...) standardGeneric("getNumStrata"))

#' @rdname gtypes.accessors
#' @export
#' 
setMethod("getNumStrata", "gtypes", function(x, ...) {
  x@data$stratum %>% 
    dplyr::n_distinct(na.rm = TRUE)
})

#' @rdname gtypes.accessors
#' @aliases getIndNames
#' @export
#' 
setGeneric("getIndNames", function(x, ...) standardGeneric("getIndNames"))

#' @rdname gtypes.accessors
#' @export
#' 
setMethod("getIndNames", "gtypes", function(x, by.strata = FALSE, ...) {
  if(by.strata) {
    x@data %>% 
      split(.$stratum) %>% 
      purrr::map(function(s) {
        s$id %>% 
          unique() %>%
          stats::na.omit() %>% 
          as.character() %>% 
          sort()
      })
  } else {
    x@data[["id"]] %>% 
      unique() %>% 
      sort()
  }
})

#' @rdname gtypes.accessors
#' @aliases getLocusNames
#' @export
#' 
setMethod("getLocusNames", "gtypes", function(x, ...) {
  x@data[["locus"]] %>% 
    unique() %>% 
    sort()
})

#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("getAlleleNames", function(x, ...) standardGeneric("getAlleleNames"))

#' @rdname gtypes.accessors
#' @aliases getAlleleNames
#' @export
#' 
setMethod("getAlleleNames", "gtypes", function(x, ...) {
  x@data %>% 
    split(.$locus) %>% 
    purrr::map(function(s) {
      s$allele %>% 
        unique() %>%
        stats::na.omit() %>% 
        sort() %>% 
        as.character()
    })
})

#' @rdname gtypes.accessors
#' @export
#' 
setGeneric("getStrataNames", function(x, ...) standardGeneric("getStrataNames"))

#' @rdname gtypes.accessors
#' @aliases getStrataNames
#' @export
#' 
setMethod("getStrataNames", "gtypes", function(x, ...) {
  x@data[["stratum"]] %>% 
    unique() %>% 
    sort()
})
          

#' #' @rdname gtypes.accessors
#' #' @export
#' #' 
#' setGeneric("ploidy", function(x, ...) standardGeneric("ploidy"))

#' @rdname gtypes.accessors
#' @aliases ploidy
#' @export
#' 
setMethod("ploidy", "gtypes", function(x, ...) x@ploidy)


#' #' @rdname gtypes.accessors
#' #' @export
#' #' 
#' setGeneric("other", function(x, ...) standardGeneric("other"))

#' @rdname gtypes.accessors
#' @aliases other
#' @export
#' 
setMethod("other", "gtypes", function(x, ...) x@other)


#' #' @rdname gtypes.accessors
#' #' @export
#' #' 
#' setGeneric("strata", function(x, ...) standardGeneric("strata"))

#' @rdname gtypes.accessors
#' @aliases strata
#' @export
#' 
setMethod("strata", "gtypes", function(x) {
  id.strata <- x@data %>% 
    dplyr::select(id, stratum) %>% 
    dplyr::distinct()
  stats::setNames(id.strata$stratum, id.strata$id)
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
  if(!is.vector(value)) {
    stop("strata must be assigned as a vector", call. = FALSE)
  }
  if(is.null(names(value))) {
    stop("strata vector must have ids for names", call. = FALSE)
  }
  missing <- !names(value) %in% x@data[["id"]]
  if(any(missing)) {
    stop(
      "the following ids are not in the gtypes object:", 
      paste(names(value)[missing], collapse = ", "),
      call. = FALSE
    )
  }
  
  x@data <- x@data %>% 
    dplyr::left_join(
      tibble::tibble(
        id = names(value),
        new = as.character(value)
      ), 
      by = "id") %>% 
    dplyr::select(id, new, locus, allele) %>% 
    dplyr::rename(stratum = new) %>% 
    data.table::as.data.table()
  
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
setGeneric("sequences", function(x, ...) standardGeneric("sequences"))

#' @rdname gtypes.accessors
#' @aliases sequences
#' @export
#' 
setMethod(
  "sequences", "gtypes", 
  function(x, as.haplotypes = TRUE, seqName = NULL, as.multidna = FALSE, ...) {
  if(is.null(x@sequences)) return(NULL)
  dna <- apex::getSequences(x@sequences, simplify = FALSE)
  if(!as.haplotypes) {
    dna <- x@data %>% 
      split(.$locus) %>% 
      purrr::map(function(l) {
        locus <- unique(l$locus)
        l <- dplyr::filter(l, !duplicated(id))
        setNames(dna[[locus]][l$allele], l$id)
      })
  }
  if(!is.null(seqName)) dna <- dna[seqName]
  if(as.multidna) as.multidna(dna) else dna
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
setMethod(
  "[", 
  signature(x = "gtypes", i = "ANY", j = "ANY", drop = "ANY"), 
  function(x, i, j, k, ..., quiet = TRUE, drop = FALSE) {
  
  # check ids (i)
  if(missing(i)) i <- TRUE
  if(is.factor(i)) i <- as.character(i)
  ids <- getIndNames(x)
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
  locs <- getLocusNames(x)
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
  st <- getStrataNames(x)
  k <- if(is.character(k)) {
    k <- unique(k)
    missing.strata <- setdiff(k, st)
    if(length(missing.strata) > 0 & !quiet) {
      missing.strata <- paste(missing.strata, collapse = ", ")
      warning("the following strata cannot be found: ", missing.strata)
    }
    intersect(k, st)
  } else st[k]

  if(length(i) == 0) stop("no samples selected")
  if(length(j) == 0) stop("no loci selected")
  if(length(k) == 0) stop("no strata selected")
  
  x@data <- x@data %>% 
    dplyr::filter(id %in% i & locus %in% j & stratum %in% k) %>% 
    data.table::as.data.table()
  if(nrow(x@data) == 0) stop("the requested indices would form an empty gtypes object")
  
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
