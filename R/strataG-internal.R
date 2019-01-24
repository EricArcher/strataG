#' @title Internal Functions
#' @description Functions intended to be used only by functions within 
#'   the strataG package.
#'   
#' @details \describe{
#'   \item{.getFileLabel}{}
#'   \item{.parseLocusNames}{}
#'   \item{.setupClusters}{}
#'   \item{.strataPairs}{}
#'   \item{.removeIdsMissingAllLoci}{}
#'   \item{.applyPerLocus}{}
#'   \item{.numericLoci}{}
#'   \item{.combineLoci}{}
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
  desc <- getDescription(g)
  label <- if(!is.null(label)) {
    label 
  } else if(!is.null(desc)) {
    desc 
  } else "strataG.gtypes"
  gsub("[[:punct:]]", ".", label)
}
  

#' @rdname strataG-internal
#' @param locus.names a vector of locus names.
#' @param ploidy integer representing the ploidy of the data.
#' @keywords internal
#' 
.expandLocusNames <- function(locus.names, ploidy) {
  if(ploidy == 1) return(locus.names)
  paste(rep(locus.names, each = ploidy), 1:ploidy, sep = ".")
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
#' @param num.cores number of cores for multithreading. 
#'   If \code{NULL}, the number used is set to the 
#'   value of \code{parallel::detectCores() - 1}.
#' @param max.cores maximum number of cores to use.
#' @keywords internal
#' 
.getNumCores <- function(num.cores, max.cores = NULL) {
  if(is.null(max.cores)) max.cores <- parallel::detectCores() - 1
  if(is.na(max.cores)) max.cores <- 1
  if(max.cores < 1) max.cores <- 1
  min(num.cores, max.cores)
}


#' @rdname strataG-internal
#' @param num.cores number of cores for multithreading. 
#'   If \code{NULL}, the number used is set to the 
#'   value of \code{parallel::detectCores() - 1}.
#' @param max.cores maximum number of cores to use.
#' @keywords internal
#' 
.setupClusters <- function(num.cores = NULL, max.cores = NULL) {
  # setup clusters
  if(is.null(num.cores)) num.cores <- parallel::detectCores() - 1
  num.cores <- .getNumCores(num.cores, max.cores)
  if(num.cores > 1) {
    cl.func <- ifelse(
      .Platform$OS.type == "windows", 
      parallel::makePSOCKcluster, 
      parallel::makeForkCluster
    )
    cl.func(num.cores)
  } else NULL
}


#' @rdname strataG-internal
#' @param g a \linkS4class{gtypes} object.
#' @keywords internal
#' 
.strataFreq <- function(g) table(getStrata(g), useNA = "no")


#' @rdname strataG-internal
#' @param g a \linkS4class{gtypes} object.
#' @keywords internal
#' 
.strataPairs <- function(g) {
  st <- getStrataNames(g)
  if(length(st) < 2) return(NULL)
  utils::combn(st, 2) %>% 
    t() %>%  
    as.data.frame(stringsAsFactors = FALSE) %>% 
    stats::setNames(c("strata.1", "strata.2"))
}


#' @rdname strataG-internal
#' @param g a \linkS4class{gtypes} object.
#' @keywords internal
#' 
.removeIdsMissingAllLoci <- function(g) {
  to.remove <- g@data %>% 
    dplyr::group_by(.data$id) %>% 
    dplyr::summarize(to.remove = all(is.na(.data$allele))) %>% 
    dplyr::filter(.data$to.remove) %>% 
    dplyr::pull(.data$id) %>% 
    as.character()

  if(length(to.remove) > 0) {
    warning(
      "The following samples are missing data for all loci and have been removed: ", 
      paste(to.remove, collapse = ", "),
      call. = FALSE
    )
    g@data <- g@data %>% 
      dplyr::filter(!.data$id %in% to.remove) %>% 
      data.table::as.data.table()
  }
  g@data <- droplevels(g@data)
  g
}


#' @rdname strataG-internal
#' @param fun a function that takes one locus column at a time.
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical - return results grouped by strata?
#' @keywords internal
#' 
.applyPerLocus <- function(fun, g, by.strata = FALSE, ...) {
  result <- if(by.strata) {
    g@data %>% 
      dplyr::group_by(.data$stratum, .data$locus) %>%  
      dplyr::summarize(value = fun(.data$allele, ...))
  } else {
    g@data %>% 
      dplyr::group_by(.data$locus) %>%  
      dplyr::summarize(value = fun(.data$allele, ...))
  }
  dplyr::ungroup(result)
}


#' @rdname strataG-internal
#' @param g a \linkS4class{gtypes} object.
#' @param min.val minimum value to start allele numbering with
#' @keywords internal
#' 
.alleles2integer <- function(g, min.val = 0) {
  g@data %>% 
    dplyr::group_by(.data$locus) %>% 
    dplyr::mutate(allele = min.val - 1 + as.integer(factor(.data$allele))) %>% 
    dplyr::ungroup() 
}


#' @rdname strataG-internal
#' @param g a \linkS4class{gtypes} object.
#' @param alleles2integer convert alleles to integers?
#' @param na.val value to replace NAs with.
#' @keywords internal
#' 
.stackedAlleles <- function(g, alleles2integer = FALSE, na.val = NULL, ...) {
  x <- if(alleles2integer) .alleles2integer(g, ...) else g@data
  if(!is.null(na.val)) {
    x <- x %>% 
      dplyr::mutate(
        allele = ifelse(is.na(.data$allele), na.val, .data$allele)
      )
  }
  x %>% 
    dplyr::arrange(.data$id, .data$locus) %>% 
    dplyr::mutate(a = rep(1:getPloidy(g), dplyr::n() / getPloidy(g))) %>% 
    tidyr::spread(.data$locus, .data$allele) %>% 
    dplyr::rename(allele = "a") %>% 
    dplyr::select(.data$id, .data$stratum, .data$allele, dplyr::everything())
}


#' @rdname strataG-internal
#' @param x a vector
#' @param sep a character for separating
#' @param sort a logical for whether to sort alleles
#' @keywords internal
#' 
.combineLoci <- function(x, sep, sort) {
  x <- as.character(x)
  if(any(is.na(x))) {
    as.character(NA)
  } else {
    x <- if(sort) sort(x) else x
    paste(x, collapse = sep)
  }
}


#' @rdname strataG-internal
#' @param g a \linkS4class{gtypes} object.
#' @keywords internal
#' 
.checkHapsLabelled <- function(g) {
  if(
    getPloidy(g) == 1 & 
    !is.null(getSequences(g)) & 
    !is.null(getOther(g, "haps.unassigned"))
  ) labelHaplotypes(g) else g
}

