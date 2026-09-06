#' @title Convert Between \code{gtypes} And \code{genind} objects.
#' @description Convert a \code{gtypes} object to a \code{genind} object 
#'   and vice-versa.
#' 
#' @param x either a \linkS4class{gtypes} or \code{genind} object
#'   to convert from.
#' @param type a character string indicating the type of marker for 
#'   \code{genind} objects: 'codom' stands for 'codominant' 
#'   (e.g. microstallites, allozymes); 'PA' stands for 'presence/absence' 
#'   markers (e.g. AFLP, RAPD). Only \code{codom} is supported; \code{PA}
#'   is rejected rather than treating marker states as codominant alleles.
#'
#' @details
#' Conversions require a single positive ploidy. Mixed-ploidy and
#' presence/absence genind objects are rejected. Original locus names are
#' retained, including dots and underscores, without normalization collisions.
#' Conversion through adegenet can discard incomplete genotypes and remove
#' samples or loci with no scored genotypes, with warnings. These conversions
#' are not lossless for partially missing calls. Allele order is not phase.
#'
#' @return A genind or gtypes object, respectively.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{initialize.gtypes}, \link{df2gtypes}, 
#'   \link{sequence2gtypes}, \link{as.data.frame.gtypes}, 
#'   \link{gtypes2loci}
#' 
#' @examples
#' data(msats.g)
#' 
#' # Convert to genind
#' gi <- gtypes2genind(msats.g)
#' gi
#' 
#' # Convert to gtypes
#' gt <- genind2gtypes(gi)
#' gt
#' 
#' @name gtypes2genind
#' @export
#' 
gtypes2genind <- function(x, type = c("codom", "PA")) {
  type <- match.arg(type)
  if(type != "codom") {
    stop("Presence/absence conversion is not supported; use codominant data.",
         call. = FALSE)
  }
  df <- as.data.frame(x, one.col = TRUE, sep = "/", strata = FALSE) |> 
    tibble::column_to_rownames("id") |> 
    as.data.frame()
  # Temporary labels avoid df2genind's locus-name normalization.
  locus.names <- stats::setNames(colnames(df), paste0("LOC", seq_len(ncol(df))))
  colnames(df) <- names(locus.names)
  
  gi <- adegenet::df2genind(
    X = df,
    sep = "/", 
    pop =  getStrata(x)[rownames(df)],
    NA.char = NA,
    ploidy = getPloidy(x),
    type = type
  )
  adegenet::locNames(gi) <- unname(locus.names[adegenet::locNames(gi)])
  adegenet::other(gi) <- getOther(x)
  gi
}


#' @rdname gtypes2genind
#' @export
#' 
genind2gtypes <- function(x) {
  if(!methods::is(x, "genind")) stop("'x' must be a genind object.", call. = FALSE)
  if(x@type != "codom") {
    stop("Presence/absence genind objects are not supported.", call. = FALSE)
  }
  ploidy <- unique(x@ploidy)
  if(length(ploidy) != 1L || is.na(ploidy) || ploidy < 1L) {
    stop("A single positive ploidy is required; mixed ploidy is not supported.",
         call. = FALSE)
  }
  locus.names <- stats::setNames(adegenet::locNames(x),
                                paste0("LOC", seq_len(adegenet::nLoc(x))))
  adegenet::locNames(x) <- names(locus.names)
  gen.mat <- adegenet::genind2df(x, usepop = TRUE, oneColPerAll = TRUE)
  gen.mat[gen.mat == "NA"] <- NA
  has.pop <- !is.null(x@pop)
  result <- df2gtypes(
    x = gen.mat,
    ploidy = ploidy,
    id.col = NULL,
    strata.col = if(has.pop) 1 else NULL,
    loc.col = if(has.pop) 2 else 1,
    schemes = x@strata,
    other = list(genind = adegenet::other(x))
  )
  # Haploid genind2df output also includes the allele-column suffix.
  old.names <- names(locus.names)
  if(ploidy == 1L) old.names <- paste0(old.names, ".1")
  lookup <- stats::setNames(unname(locus.names), old.names)
  data.table::set(result@data, j = "locus",
                  value = unname(lookup[result@data$locus]))
  data.table::setkeyv(result@data, c("id", "stratum", "locus"))
  result
}
