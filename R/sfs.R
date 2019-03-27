#' @title Site Frequency Spectrum
#' @description Calculate the SFS from a data frame of SNP genotypes
#'   
#' @param x a data frame of SNP genotypes where the first two columns are id and
#'   strata designations and SNPs start on the third column. SNP genotypes are
#'   coded as integers where 0 and 2 are the major and minor homozygotes and 1
#'   is the heterozygote.
#' @param joint return the joint SFS if more than 1 strata is available?
#' @param ambig.maf treat the MAF designation as potentially ambiguous? If
#'   `TRUE`, then loci with a global MAF = 0.5 will contribute 0.5 to the 
#'   SFS rather than 1. If `FALSE`, then the MAF is assumed to be unambiguous
#'   for all loci regardless of the observed frequency.
#' @param fsc.dimnames format matrix dimnames for fastsimcoal2? If `TRUE`, then
#'   row and column names will be prefixed with the deme number (e.g., "d0_") 
#'   that they represent.
#'    
#' @return Either a named vector of the 1D site frequency 
#'   spectrum if \code{joint = FALSE}, or a list of matrices of the 2D
#'   site frequency spectra for pairs of strata if \code{joint = TRUE}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
sfs <- function(x, joint = FALSE, ambig.maf = TRUE, fsc.dimnames = FALSE) {
  if(!is.data.frame(x)) stop("'x' must be a data frame.")
  st <- x[, 2]
  x <- x[, -(1:2), drop = FALSE]
  global.count <- colSums(x)
  wt <- if(ambig.maf) ifelse(global.count == nrow(x), 0.5, 1) else NULL
  if(!joint | dplyr::n_distinct(st) == 1) {
    loc.freq <- if(ambig.maf) {
      tapply(wt, global.count, sum)
    } else {
      table(global.count)
    }
    sfs.1d <- stats::setNames(vector("numeric", nrow(x) + 1), 0:nrow(x))
    sfs.1d[names(loc.freq)] <- loc.freq
    if(fsc.dimnames) names(sfs.1d) <- paste0("d0_", names(sfs.1d))
    sfs.1d
  } else {
    counts <- lapply(split(x, st), colSums)
    n <- table(st)
    combn(unique(st), 2, function(pair) {
      i <- pair[1]
      j <- pair[2]
      freqs <- if(ambig.maf) {
        tapply(wt, list(counts[[i]], counts[[j]]), sum) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column("freq1") %>% 
          tidyr::gather("freq2", "count", -.data$freq1) %>% 
          dplyr::mutate(count = ifelse(is.na(.data$count), 0, .data$count))
      } else {
        as.data.frame(table(counts[[i]], counts[[j]]), stringsAsFactors = FALSE)
      }
      nrow <- (2 * n[i]) + 1
      ncol <- (2 * n[j]) + 1
      sfs.2d <- matrix(0, nrow = nrow, ncol = ncol)
      dimnames(sfs.2d) <- stats::setNames(list(0:(nrow - 1), 0:(ncol - 1)), pair)
      for(r in 1:nrow(freqs)) sfs.2d[freqs[r, 1], freqs[r, 2]] <- freqs[r, 3]
      
      
      
      t(sfs.2d)
    }, simplify = FALSE)
  }
}