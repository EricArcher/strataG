#' @title \code{gtypes} Constructor
#' @description Create a new \linkS4class{gtypes} object using 
#'   \code{new("gtypes", ...)}, where '\code{...}' are arguments 
#'   documented below.
#' 
#' @param .Object the object skeleton, automatically generated when 
#'   calling \code{new}.
#' @param gen.data a vector, matrix, or data.frame containing the alleles 
#'   at each locus. See below for more details.
#' @param ploidy ploidy of the loci.
#' @param ind.names an optional vector of individual sample names.
#' @param sequences an optional \linkS4class{multidna} 
#'   object containing sequences represented by each locus.
#' @param strata an optional stratification scheme from \code{schemes}.
#' @param schemes an optional data.frame of stratification schemes.
#' @param description an optional description for the object.
#' @param other other optional information to include.
#' @param remove.sequences logical. If \code{TRUE} any sequences not referenced 
#'
#' @details
#' For multi-allele loci, the \code{gen.data} argument should be 
#' formatted such that every consecutive \code{ploidy} columns represent 
#' alleles of one locus. Locus names are taken from the column names in 
#' \code{gen.data} and should be formatted with the same root locus name, with 
#' unique suffixes representing alleles (e.g., for Locus1234: Locus1234.1 
#' and Locus1234.2, or Locus1234_A and Locus1234_B). \cr
#' If \code{gen.data} is a vector it is assumed to represent haplotypes of a 
#' haploid marker.
#' Sample names can be either in the rownames of \code{gen.data} or given 
#' separately in \code{ind.names}. If \code{ind.names} are provided, these are 
#' used in lieu of rownames in \code{gen.data}. 
#' If \code{schemes} has a column named '\code{id}', it will be used to match
#' to sample names in \code{gen.data}. Otherwise, if rownames are present in
#' \code{schemes}, a column named '\code{id}' will be created from them.
#' If sequences are provided in \code{sequences}, then they should be named 
#' and match values in the haplotype column in \code{gen.data}. If multiple 
#' genes are given as a \linkS4class{multidna} object, it is assumed that they
#' are in the same order as the columns in \code{gen.data}.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @seealso \link{df2gtypes}, \link{sequence2gtypes}
#  \link{gtypes2genind}, \link{gtypes2loci}
#'
#' @aliases gtypes.initialize initialize.gtypes new
#' 
methods::setMethod(
  "initialize", 
  "gtypes",
  function(.Object, gen.data, ploidy, ind.names = NULL,
           sequences = NULL, strata = NULL, schemes = NULL,
           description = NULL, other = NULL,
           remove.sequences = FALSE) {
  
  if(is.null(gen.data) | is.null(ploidy)) return(.Object)
            
  # check gen.data
  if(is.vector(gen.data)) {
    gen.data <- cbind(gen.data)
    colnames(gen.data) <- "Haplotype"
  }
  if(!(is.matrix(gen.data) | is.data.frame(gen.data))) {
    stop("'gen.data' is not a vector, matrix, or data.frame", call. = FALSE)
  }
  gen.data <- as.matrix(gen.data)
  
  # check ploidy
  ploidy <- as.integer(ploidy)
  if(ncol(gen.data) %% ploidy != 0) {
    stop(
      "the number of columns in 'gen.data' is not a multiple of 'ploidy'",
      call. = FALSE
    )
  }
  if(ploidy > 1 & !is.null(sequences)) {
    stop(
      "'sequences' can't be present if 'ploidy' is greater than 1", 
      call. = FALSE
    )
  }
  
  # check ind.names
  ind.names <- if(!is.null(ind.names)) {
    if(!is.vector(ind.names)) stop("'ind.names' must be a vector")
    if(length(ind.names) != nrow(gen.data)) {
      stop(
        "the length of 'ind.names' must equal the number of rows in 'gen.data'",
        call. = FALSE
      )
    }
    as.character(ind.names)
  } else {
    if(!is.null(rownames(gen.data))) rownames(gen.data) else 1:nrow(gen.data)
  }
  if(any(duplicated(ind.names))) {
    dup.names <- unique(ind.names[duplicated(ind.names)])
    dup.names <- paste(dup.names, collapse = ", ")
    stop("there are duplicated individual names: ", dup.names, call. = FALSE)
  }
  rownames(gen.data) <- ind.names
  
  # check strata
  if(!is.null(strata)) {
    if(is.null(names(strata))) {
      if(length(strata) == length(ind.names)) names(strata) <- ind.names
    } 
  } else strata <- rep("Default", nrow(gen.data))
  if(length(strata) != nrow(gen.data)) {
    warning("the length of 'strata' is not the same as the number of individuals. strata will be recycled.")
  }
  
  # check schemes
  if(!is.null(schemes)) {
    # check that schemes is a data.frame
    if(!(is.matrix(schemes) | is.data.frame(schemes))) {
      stop("'schemes' is not a matrix or data.frame", call. = FALSE)
    } 
    schemes <- as.data.frame(schemes)

    # check that 'id' column in schemes exists
    if(!"id" %in% colnames(schemes)) {
      if(is.null(rownames(schemes))) {
        if(nrow(schemes) != nrow(gen.data)) {
          stop(
            "'schemes' doesn't have an 'id' column or rownames and is not as long as 'gen.data'",
            call. = FALSE
          )
        } else {
          schemes$id <- 1:nrow(schemes)
        }
      } else {
        schemes$id <- rownames(schemes)
      }
    }
    rownames(schemes) <- NULL
    schemes <- dplyr::select(schemes, .data$id, dplyr::everything())

    # check that ids in schemes can be found
    if(length(intersect(schemes$id, rownames(gen.data))) == 0) {
      stop(
        "no ids in 'schemes' are in 'gen.data' or 'ind.names'", 
        call. = FALSE
      )
    }
    
    schemes$id <- as.character(schemes$id)
  }
  
  # check description
  if(is.null(description)) {
    description <- paste(
      "gtypes created on", format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
  }
  
  # create locus names
  nloc <- ncol(gen.data) / ploidy
  if(is.null(colnames(gen.data))) {
    # return generic names if no colnames assigned
    nums <- formatC(1:nloc, digits = floor(log10(nloc)), flag = "0")
    generic.locus.names <- paste0("Locus", "_", nums)
    colnames(gen.data) <- .expandLocusNames(generic.locus.names, ploidy)
  } 
  locus.names.lookup <- stats::setNames(
    rep(.parseLocusNames(colnames(gen.data), ploidy), each = ploidy),
    colnames(gen.data)
  )
  
  # check sequences
  if(!is.null(sequences)) {
    sequences <- as.multidna(sequences)
    if(getNumLoci(sequences) != ncol(gen.data)) {
      stop(
        "the number of genes in 'sequences' is not equal to the number of loci",
        call. = FALSE
      )
    }
    setLocusNames(sequences) <- colnames(gen.data)
    for(loc in colnames(gen.data)) {
      haps <- stats::na.omit(unique(as.character(gen.data[, loc])))
      seq.names <- apex::getSequenceNames(sequences)[[loc]]
      missing <- setdiff(haps, seq.names)
      if(length(missing) > 0) {
        stop(
          "the following haplotypes can't be found in sequences for locus '", 
          loc, "': ", paste(missing, collapse = ", "), 
          call. = FALSE
        )
      }
    }
  }
  
  gen.data <- cbind(
    id = rownames(gen.data), 
    stratum = as.character(strata), 
    gen.data
  ) %>% 
    as.data.frame(stringsAsFactors = FALSE) %>% 
    tidyr::gather("locus", "allele", -.data$id, -.data$stratum) %>% 
    dplyr::mutate(locus = locus.names.lookup[as.character(.data$locus)])
  
  data.table::setDT(gen.data, key = c("id", "stratum", "locus"))
  
  # create and return gtypes object
  g <- .Object
  g@data <- gen.data
  g@ploidy <- ploidy
  g@sequences <- sequences
  g@schemes <- schemes
  g@description <- description
  g@other <- if(is.null(other)) list() else other

  # Check for samples missing data for all loci
  g <- .removeIdsMissingAllLoci(g)

  # Remove unreferenced sequences
  if(remove.sequences) g <- removeSequences(g)

  g
})
