#' @title Site Frequency Spectrum
#' @description Calculate the SFS from a data frame of SNP genotypes
#'   
#' @param x a data frame of SNP genotypes where the first two columns are id and
#'   strata designations and SNPs start on the third column. SNP genotypes are
#'   coded as integers where 0 and 2 are the major and minor homozygotes and 1
#'   is the heterozygote.
#' @param fsc.dimnames format matrix dimnames for fastsimcoal2? If \code{TRUE},
#'   then row and column names will be prefixed with the deme number (e.g.,
#'   "d0_") that they represent.
#' @param strata.col column number that strata designations are in.
#' @param locus.col column number that loci start in. All columns after this are
#'   assumed to be loci.
#' @param sort.strata if \code{joint = TRUE}, are strata to be sorted
#'   alphabetically? If \code{FALSE} then strata are taken in the order found in
#'   \code{strata.col}.
#' @param na.action action to take if genotypes are missing for some samples. 
#'   If \code{"fail"}, an error is thrown if any genotypes are missing. If 
#'   \code{"filter"}, SNPs with missing genotypes are removed.
#'    
#' @return A list of the marginal (1D) and joint (2D) site frequency spectra. 
#'   If only one stratum is present, then \code{$marginal} will be  \code{NULL}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
sfs <- function(x, strata.col = 2, locus.col = 3, fsc.dimnames = TRUE, 
                sort.strata = TRUE, na.action = c("fail", "filter")) {
  if(!is.data.frame(x)) stop("'x' must be a data frame.")
  if(strata.col >= locus.col) stop("'strata.col' can't be >= 'locus.col'.")
  
  st <- x[, strata.col]
  st.n <- table(st)    
  st.order <- unique(st)
  if(sort.strata) st.order <- st.order[order(nchar(st.order), st.order)]
  
  x <- x[, locus.col:ncol(x), drop = FALSE]
  has.missing <- sapply(x, function(loc) any(is.na(loc)))
  if(any(has.missing)) {
    if(match.arg(na.action) == "fail") {
      stop("'x' has missing data. Can't compute SFS with na.action = 'fail'.")
    } 
    to.keep <- apply(x, 2, function(loc) all(!is.na(loc)))
    if(sum(to.keep) == 0) {
      stop("Can't compute SFS because all loci have missing data.")
    } 
    warning(
      ncol(x) - sum(to.keep), " loci had missing data. ",
      "Retained ", sum(to.keep), " loci.", 
      call. = FALSE
    )
    x <- x[, to.keep]
  }
  
  st.x <- split(x, st)
  st.count <- sapply(st.x, colSums, simplify = FALSE)
  is.ambig <- colSums(x) == nrow(x)
  
  marginal <- sapply(st.order, function(k) {
    sfs.vec <- vector("numeric", length = (st.n[k] * 2) + 1)
    for(loc in 1:length(st.count[[k]])) {
      i1 <- st.count[[k]][loc] + 1
      if(is.ambig[loc]) {
        i2 <- length(sfs.vec) - i1 + 1
        sfs.vec[i1] <- sfs.vec[i1] + 0.5
        if(i1 != i2) sfs.vec[i2] <- sfs.vec[i2] + 0.5
      } else {
        sfs.vec[i1] <- sfs.vec[i1] + 1
      }
    }
    names(sfs.vec) <- 0:(length(sfs.vec) - 1)
    if(fsc.dimnames) {
      names(sfs.vec) <- paste0("d", match(k, st.order) - 1, "_", names(sfs.vec))
    }
    sfs.vec
  }, simplify = FALSE)
  
  joint <- if(length(st.order) == 1) NULL else {
    utils::combn(st.order, 2, function(pair) {
      k1 <- pair[1]
      k2 <- pair[2]
      count1 <- st.count[[k1]]
      count2 <- st.count[[k2]]
      sfs.mat <- matrix(0, nrow = (2 * st.n[k1]) + 1, ncol = (2 * st.n[k2]) + 1)
      for(loc in 1:length(count1)) {
        i1 <- count1[loc] + 1
        j1 <- count2[loc] + 1
        if(is.ambig[loc]) {
          i2 <- nrow(sfs.mat) - i1 + 1
          j2 <- ncol(sfs.mat) - j1 + 1
          sfs.mat[i1, j1] <- sfs.mat[i1, j1] + 0.5
          sfs.mat[i2, j2] <- sfs.mat[i2, j2] + 0.5
        } else {
         sfs.mat[i1, j1] <- sfs.mat[i1, j1] + 1
        }
      }
      rownames(sfs.mat) <- 0:(nrow(sfs.mat) - 1)
      colnames(sfs.mat) <- 0:(ncol(sfs.mat) - 1)
      if(fsc.dimnames) {
        row.deme <- paste0("d",  match(k1, st.order) - 1, "_")
        col.deme <- paste0("d",  match(k2, st.order) - 1, "_")
        rownames(sfs.mat) <- paste0(row.deme, rownames(sfs.mat))
        colnames(sfs.mat) <- paste0(col.deme, colnames(sfs.mat))
      }
      names(dimnames(sfs.mat)) <- pair
      t(sfs.mat)
    }, simplify = FALSE)
  }
  
  list(marginal = marginal, joint = joint)
}