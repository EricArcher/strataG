#' @title Convert Between \code{gtypes} And \code{loci} objects.
#' @description Convert a \code{gtypes} object to a \code{\link[pegas]{loci}} object.
#' 
#' @param x a \linkS4class{gtypes} or \code{loci} formatted object.
#' @param sep a single literal character separating alleles. Forward conversion
#'   returns canonical slash-separated pegas genotypes. Reverse conversion
#'   accepts canonical slash/pipe separators or the supplied separator.
#' @param description a label for the \code{gtypes} object (optional).
#'  
#' @details
#' Only diploid data are supported. Missing genotypes remain actual NA values;
#' zero is not used as a missing-data sentinel. Partially missing calls become
#' missing genotypes during forward conversion. Missing population labels do
#' not cause samples to be removed. Reverse conversion uses the population
#' column, wherever it occurs, or the default stratum when it is absent.
#' The gtypes constructor may remove samples missing every genotype, with a
#' warning. Haploid and other non-diploid calls are rejected; ploidy cannot be
#' inferred from an entirely missing locus.
#'
#' @return A loci or gtypes object, respectively.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{initialize.gtypes}, \link{df2gtypes}, 
#'   \link{sequence2gtypes}, \link{as.data.frame.gtypes}, 
#'   \link{gtypes2genind}
#' 
#' @examples
#' data(msats.g)
#' 
#' # Convert to loci
#' lc <- gtypes2loci(msats.g)  
#' lc  
#' 
#' # Convert to gtypes
#' gt <- loci2gtypes(lc)
#' gt 
#' 
#' @name gtypes2loci 
#' @export
#' 
gtypes2loci <- function(x, sep = "/") {
  if(!is.gtypes(x)) stop("'x' must be a gtypes object")
  if(getPloidy(x) != 2L) stop("Only diploid data are supported.", call. = FALSE)
  if(!is.character(sep) || length(sep) != 1L || is.na(sep) || nchar(sep) != 1L) {
    stop("'sep' must be one non-missing character.", call. = FALSE)
  }
  alleles <- x@data$allele
  if(any(vapply(unique(c(sep, "/", "|")), function(s) {
    any(grepl(s, alleles[!is.na(alleles)], fixed = TRUE))
  }, logical(1)))) {
    stop("Allele labels cannot contain genotype separators.", call. = FALSE)
  }
  df <- as.data.frame(x, one.col = TRUE, sep = sep)
  rownames(df) <- df$id
  df$id <- NULL
  # pegas loci use canonical slash-separated genotypes.
  for(i in seq.int(2L, ncol(df))) {
    df[[i]] <- factor(gsub(sep, "/", df[[i]], fixed = TRUE))
  }
  pegas::as.loci(df, allele.sep = "/|", col.pop = 1)
}

#' @rdname gtypes2loci
#' @export
#'
loci2gtypes <- function(x, description = NULL, sep = "/") {
  if(!inherits(x, "loci")) stop("'x' must be a loci object")
  if(!is.character(sep) || length(sep) != 1L || is.na(sep) || nchar(sep) != 1L) {
    stop("'sep' must be one non-missing character.", call. = FALSE)
  }
  locicol <- attr(x, "locicol")
  if(!is.numeric(locicol) || !length(locicol) || anyNA(locicol) ||
     any(locicol != as.integer(locicol)) || anyDuplicated(locicol) ||
     any(locicol < 1L | locicol > ncol(x))) {
    stop("The loci object must identify valid genotype columns.", call. = FALSE)
  }
  df <- as.data.frame(x)
  popcol <- which(names(df) == "population" & !seq_along(df) %in% locicol)
  if(length(popcol) > 1L) stop("Multiple population columns found.", call. = FALSE)
  pop <- if(length(popcol)) as.character(df[[popcol]]) else rep("Default", nrow(df))
  columns <- lapply(locicol, function(i) {
    values <- as.character(df[[i]])
    result <- matrix(NA_character_, nrow(df), 2L)
    for(j in which(!is.na(values))) {
      value <- values[j]
      # as.loci normalizes separators; also accept explicitly separated input.
      separator <- if(grepl("/", value, fixed = TRUE)) "/" else
        if(grepl("|", value, fixed = TRUE)) "|" else sep
      parts <- strsplit(value, separator, fixed = TRUE)[[1]]
      if(length(parts) != 2L || any(!nzchar(parts)) ||
         startsWith(value, separator) || endsWith(value, separator)) {
        stop("Only diploid genotypes with two separated alleles are supported.",
             call. = FALSE)
      }
      result[j, ] <- parts
    }
    result
  })
  mat <- do.call(cbind, columns)
  # Temporary names avoid data.frame name repair and locus-prefix parsing.
  keys <- paste0("LOC", seq_along(locicol))
  colnames(mat) <- paste0(rep(keys, each = 2L), ".", rep(1:2, length(keys)))
  input <- data.frame(id = rownames(df), pop = pop, mat, check.names = FALSE)
  result <- df2gtypes(input, ploidy = 2, description = description)
  lookup <- stats::setNames(names(df)[locicol], keys)
  data.table::set(result@data, j = "locus",
                  value = unname(lookup[result@data$locus]))
  data.table::setkeyv(result@data, c("id", "stratum", "locus"))
  result
}
