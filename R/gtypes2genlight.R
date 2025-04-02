#' @title Convert Between \code{gtypes} And \code{genlight} objects.
#' @description Convert a \code{gtypes} object to a \code{genlight} object 
#'   and vice-versa.
#' 
#' @param x either a \linkS4class{gtypes} or \code{genlight} object
#'   to convert from.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{initialize.gtypes}, \link{df2gtypes}, 
#'   \link{sequence2gtypes}, \link{as.data.frame.gtypes}, 
#'   \link{gtypes2loci}, \link{gtypes2genind}
#' 
#' @examples
#' data(msats.g)
#' 
#' # Create simple simulated SNPs
#' gl1 <- adegenet::glSim(n.ind = 100, n.snp.nonstruc = 1000, ploidy = 2)
#' gl1
#' 
#' # Convert to gtypes
#' gt <- genlight2gtypes(gl1)
#' gt
#' 
#' # Convert back to genlight
#' gl2 <- gtypes2genlight(gt)
#' gl2
#' 
#' @name gtypes2genlight
#' @export
#' 
gtypes2genlight <- function(x) {
  gl <- adegenet::as.genlight(
    as.data.frame(x, ids = FALSE, strata = FALSE, coded.snps = TRUE)
  )
  adegenet::locNames(gl) <- getLociNames(x)
  adegenet::indNames(gl) <- getIndNames(x)
  adegenet::pop(gl) <- getStrata(x)
  adegenet::other(gl) <- getOther(x)
  gl
}

#' @rdname gtypes2genlight
#' @export
#' 
genlight2gtypes <- function(x) {
  if(!inherits(x, "genlight")) stop("'x' must be a genlight object")
  genotypes <- list(c("A", "A"), c("A", "G"), c("G", "G"), as.character(c(NA, NA)))
  gen.mat <- do.call(
    cbind, 
    lapply(as.data.frame(x), function(num.alt) {
      num.alt <- ifelse(is.na(num.alt), 4, num.alt + 1)
      do.call(rbind, genotypes[num.alt])
    })
  )
  loci <- x@loc.names
  if(is.null(loci)) loci <- paste0("L", 1:adegenet::nLoc(x))
  colnames(gen.mat) <- paste0(rep(loci, each = 2), ".", 1:2)
  has.pop <- !is.null(x@pop)
  gen.mat <- cbind(
    id = 1:nrow(gen.mat),
    strata = if(has.pop) x@pop else "Default",
    as.data.frame(gen.mat)
  )
  df2gtypes(
    x = gen.mat,
    ploidy = 2,
    schemes = x@strata,
    other = list(genind = adegenet::other(x))
  )
}
