#' @title \code{gtypes} Accessors
#' @description Accessors for slots in \linkS4class{gtypes} objects.
#'
#' @param x a \linkS4class{gtypes} object.
#' @param by.strata logical - return results by strata?
#' @param seqName the name (or number) of a set of sequences from the
#'   \code{@@sequences} slot to return.
#' @param as.haplotypes return sequences as haplotypes? If \code{TRUE}, contents
#'   of \code{@@sequences} slot are returned. If \code{FALSE}, one sequence per
#'   individual is returned.
#' @param as.multidna return sequences as a \linkS4class{multidna} object? If
#'   \code{FALSE}, sequences are returned as a list.
#' @param i,j,k subsetting slots for individuals (\code{i}), loci (\code{j}), or
#'   strata (\code{k}). See Details for more information.
#' @param quiet suppress warnings about unmatched requested individuals, loci,
#'   or strata?
#' @param drop if \code{TRUE} the return object will have unused sequences
#'   removed.
#' @param ... other arguments passed from generics (ignored).
#' @param value value being assigned by accessor.
#' @param name name of the value going into the \code{other} list.
#'
#' @details Indexing a \code{gtypes} object with integers, characters, or
#' logicals with the \code{[} operator follows the same rules as normal indexing
#' in R. The order that individuals, loci, and strata are chosen is the order
#' returned by \code{getIndNames}, \code{getLocNames}, and \code{getStrataNames}
#' respectively. If unstratified samples are present, they can be selected as a
#' group either by including \code{NA} in the character or numeric vector of the
#' \code{k} slot, or by providing a logical vector based on
#' \code{is.na(strata(g))} to the \code{i} slot.
#'
#' @return \describe{ \item{nInd}{number of individuals} \item{nLoc}{number of
#' loci} \item{nStrata}{number of strata} \item{indNames}{vector of
#' individual/sample names} \item{locNames}{vector of locus names}
#' \item{strataNames}{vector of strata names for current scheme}
#' \item{ploidy}{number of alleles at each locus} \item{other}{contents of
#' \code{@@other} slot} \item{strata}{return or modify the current
#' stratification} \item{schemes}{return or modify the current stratification
#' schemes} \item{alleleNames}{return a list of alleles at each locus}
#' \item{sequences}{return the \linkS4class{multidna} object in the
#' \code{@@sequences} slot. See \code{\link[apex]{getSequences}} to extract
#' individual genes or sequences from this object} \item{description}{return the
#' object's description} }
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' #--- create a diploid (microsatellite) gtypes object
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#'
#' getNumStrata(msats.g)
#' getStrataNames(msats.g)
#' getNumLoci(msats.g)
#' getLociNames(msats.g)
#'
#' # reassign all samples to two randomly chosen strata
#' new.strata <- sample(c("A", "B"), getNumInd(msats.g), rep = TRUE)
#' names(new.strata) <- getIndNames(msats.g)
#' setStrata(msats.g) <- new.strata
#' msats.g
#'
#'
#' #--- a sequence example
#' library(ape)
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' wood.g <- sequence2gtypes(x)
#' new.strata <- sample(c("A", "B"), getNumInd(wood.g), rep = TRUE)
#' names(new.strata) <- getIndNames(wood.g)
#' setStrata(wood.g) <- new.strata
#' wood.g
#'
#' # get the multidna sequence object
#' multi.seqs <- getSequences(wood.g, as.multidna = TRUE)
#' class(multi.seqs) # "multidna"
#'
#' # get a list of DNAbin objects
#' dnabin.list <- getSequences(wood.g)
#' class(dnabin.list) # "list"
#'
#' # get a DNAbin object of the first locus
#' dnabin.1 <- getSequences(wood.g)[[1]]
#' class(dnabin.1) # "DNAbin"
#'
#' # getting and setting values in the `other` slot:
#' getOther(dloop.g)
#'
#' setOther(dloop.g, "timestamp") <- timestamp()
#' setOther(dloop.g, "Author") <- "Hoban Washburne"
#'
#' getOther(dloop.g)
#' getOther(dloop.g, "timestamp")
#'
#' setOther(dloop.g, "Author") <- NULL
#' getOther(dloop.g)
#'
#'
#' @name gtypes.accessors
#' @aliases accessors
#'   
NULL

#' @rdname gtypes.accessors
#' @aliases getNumInd
#' @export
#' 
methods::setMethod("getNumInd", "gtypes", function(x, by.strata = FALSE, ...) {
  if(by.strata) {
    x@data %>% 
      dplyr::group_by(.data$stratum) %>% 
      dplyr::summarize(num.ind = dplyr::n_distinct(.data$id)) %>% 
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
methods::setMethod("getNumLoci", "gtypes", function(x, ...) {
  dplyr::n_distinct(x@data$locus, na.rm = TRUE)
})

#' @rdname gtypes.accessors
#' @aliases getNumStrata
#' @export
#' 
methods::setGeneric("getNumStrata", function(x, ...) standardGeneric("getNumStrata"))

#' @rdname gtypes.accessors
#' @export
#' 
methods::setMethod("getNumStrata", "gtypes", function(x, ...) {
  dplyr::n_distinct(x@data$stratum, na.rm = TRUE)
})

#' @rdname gtypes.accessors
#' @aliases getIndNames
#' @export
#' 
methods::setGeneric("getIndNames", function(x, ...) standardGeneric("getIndNames"))

#' @rdname gtypes.accessors
#' @export
#' 
methods::setMethod("getIndNames", "gtypes", function(x, by.strata = FALSE, ...) {
  if(by.strata) {
    split(x@data$stratum) %>% 
      purrr::map(function(s) {
        s$id %>% 
          unique() %>%
          stats::na.omit() %>% 
          as.character() %>% 
          sort()
      })
  } else {
    sort(unique(x@data[["id"]]))
  }
})

#' @rdname gtypes.accessors
#' @export
#' 
methods::setGeneric("getLociNames", function(x, ...) standardGeneric("getLociNames"))

#' @rdname gtypes.accessors
#' @aliases getLociNames
#' @export
#' 
methods::setMethod("getLociNames", "gtypes", function(x, ...) {
  sort(unique(x@data[["locus"]]))
})

#' @rdname gtypes.accessors
#' @export
#' 
methods::setGeneric("getAlleleNames", function(x, ...) standardGeneric("getAlleleNames"))

#' @rdname gtypes.accessors
#' @aliases getAlleleNames
#' @export
#' 
methods::setMethod("getAlleleNames", "gtypes", function(x, ...) {
  split(x@data, x@data$locus) %>% 
    purrr::map(function(s) {
      as.character(sort(stats::na.omit(unique(s$allele))))
    })
})

#' @rdname gtypes.accessors
#' @export
#' 
methods::setGeneric("getStrataNames", function(x, ...) standardGeneric("getStrataNames"))

#' @rdname gtypes.accessors
#' @aliases getStrataNames
#' @export
#' 
methods::setMethod("getStrataNames", "gtypes", function(x, ...) {
  sort(unique(x@data[["stratum"]]))
})
          

#' @rdname gtypes.accessors
#' @export
#'
methods::setGeneric("getPloidy", function(x, ...) standardGeneric("getPloidy"))

#' @rdname gtypes.accessors
#' @aliases getPloidy
#' @export
#' 
methods::setMethod("getPloidy", "gtypes", function(x, ...) x@ploidy)


#' @rdname gtypes.accessors
#' @export
#'
methods::setGeneric("getStrata", function(x, ...) standardGeneric("getStrata"))

#' @rdname gtypes.accessors
#' @aliases getStrata
#' @export
#' 
methods::setMethod("getStrata", "gtypes", function(x) {
  id.strata <- x@data %>% 
    dplyr::select(.data$id, .data$stratum) %>% 
    dplyr::distinct()
  stats::setNames(id.strata$stratum, id.strata$id)
})


#' @rdname gtypes.accessors
#' @export
#' 
methods::setGeneric("setStrata<-", function(x, value) standardGeneric("setStrata<-"))

#' @rdname gtypes.accessors
#' @aliases setStrata
#' @export
#' 
methods::setMethod("setStrata<-", "gtypes", function(x, value) {
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
      tibble::tibble(id = names(value), .new = as.character(value)), 
      by = "id"
    ) %>% 
    dplyr::select(.data$id, .data$.new, .data$locus, .data$allele) %>% 
    dplyr::rename(stratum = .data$.new) %>% 
    data.table::as.data.table()
  
  methods::validObject(x)
  x
})


#' @rdname gtypes.accessors
#' @export
#' 
methods::setGeneric("getSchemes", function(x, ...) standardGeneric("getSchemes"))

#' @rdname gtypes.accessors
#' @aliases getSchemes
#' @export
#' 
methods::setMethod("getSchemes", "gtypes", function(x, ...) x@schemes)

#' @rdname gtypes.accessors
#' @export
#' 
methods::setGeneric("setSchemes<-", function(x, value) standardGeneric("setSchemes<-"))

#' @rdname gtypes.accessors
#' @aliases setSchemes
#' @export
#' 
methods::setMethod("setSchemes<-", "gtypes", function(x, value) {
  x@schemes <- value
  methods::validObject(x)
  x
})



#' @rdname gtypes.accessors
#' @export
#' 
methods::setGeneric("getSequences", function(x, ...) standardGeneric("getSequences"))

#' @rdname gtypes.accessors
#' @aliases getSequences
#' @export
#' 
methods::setMethod(
  "getSequences", 
  "gtypes", 
  function(x, as.haplotypes = TRUE, seqName = NULL, as.multidna = FALSE, ...) {
    if(is.null(x@sequences)) return(NULL)
    dna <- apex::getSequences(x@sequences, simplify = FALSE)
    if(!as.haplotypes) {
      dna <- purrr::map(
        split(x@data, x@data$locus),
        function(l) {
          locus <- unique(l$locus)
          l <- dplyr::filter(l, !duplicated(.data$id))
          stats::setNames(dna[[locus]][l$allele], l$id)
        }
      )
    }
    if(!is.null(seqName)) dna <- dna[seqName]
    if(as.multidna) as.multidna(dna) else dna
})


#' @rdname gtypes.accessors
#' @export
#' 
methods::setGeneric("getDescription", function(x, ...) standardGeneric("getDescription"))

#' @rdname gtypes.accessors
#' @aliases getDescription
#' @export
#' 
methods::setMethod("getDescription", "gtypes", function(x, ...) x@description)

#' @rdname gtypes.accessors
#' @export
#' 
methods::setGeneric("setDescription<-", function(x, value) standardGeneric("setDescription<-"))

#' @rdname gtypes.accessors
#' @aliases setDescription
#' @export
#' 
methods::setMethod("setDescription<-", "gtypes", function(x, value) {
  x@description <- value
  methods::validObject(x)
  x
})



#' @rdname gtypes.accessors
#' @export
#'
methods::setGeneric("getOther", function(x, ...) standardGeneric("getOther"))

#' @rdname gtypes.accessors
#' @aliases getOther
#' @export
#' 
methods::setMethod("getOther", "gtypes", function(x, value = NULL, ...) {
  if(is.null(value)) {
    x@other 
  } else if(length(value) == 1) {
    x@other[[value]]
  } else {
    x@other[value]
  }
})

#' @rdname gtypes.accessors
#' @export
#' 
methods::setGeneric("setOther<-", function(x, name, value) standardGeneric("setOther<-"))

#' @rdname gtypes.accessors
#' @aliases setOther
#' @export
#' 
methods::setMethod("setOther<-", "gtypes", function(x, name, value) {
  x@other[[name]] <- value
  methods::validObject(x)
  x
})



#' @rdname gtypes.accessors
#' @aliases index subset
#' @export
#' 
methods::setMethod(
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
  locs <- getLociNames(x)
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

  # subset data
  if(length(i) == 0) stop("no samples selected")
  if(length(j) == 0) stop("no loci selected")
  if(length(k) == 0) stop("no strata selected")

  x@data <- x@data[id %in% i & locus %in% j & stratum %in% k]
  if(nrow(x@data) == 0) {
    stop("the requested indices would form an empty gtypes object")
  }
  
  # filter sequences
  if(!is.null(x@sequences)) {
    j.seqs <- getSequences(
      x@sequences, loci = j, simplify = FALSE, 
      exclude.gap.only = FALSE
    )
    x@sequences <- methods::new("multidna", j.seqs)
  }
  
  # Check for samples missing data for all loci
  x <- .removeIdsMissingAllLoci(x)
  
  if(drop) x <- removeSequences(x)
  return(x)
})
