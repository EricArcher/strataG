#' @title Site Frequency Spectrum
#' @description Calculate the SFS from a \linkS4class{gtypes} object for 
#'   all individuals or pairs of strata.
#'   
#' @param g a \linkS4class{gtypes} object.
#' @param ref.allele an optional vector of reference alleles for each SNP. 
#'   See \code{\link{as.data.frame.gtypes}} for more information.
#' @param joint a logical specifying whether or not to compute joint (2D) 
#'   SFS.
#' @param strata a vector of strata to use to compute joint SFS. 
#'   If `NULL` (default), then all pairs of strata in `g` are used.
#' @param simplify a logical specifying whether to return a matrix if only 
#'   a single joint SFS is computed (one pair of strata).
#' @param na.action action to take if genotypes are missing for some samples. 
#'   If \code{"fail"}, an error is thrown if any genotypes are missing. If 
#'   \code{"filter"}, SNPs with missing genotypes are removed. If 
#'   \code{"resample"}, the minimum number of genotyped individuals across 
#'   all SNPs is first determined, then this number of non-missing genotypes 
#'   are randomly selected for each SNP to calculate the SFS.
#' 
#' @return Either a named vector of the 1D site frequency 
#'   spectrum if \code{joint = FALSE}, or a list of matrices of the 2D
#'   site frequency spectra for pairs of strata if \code{joint = TRUE}.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
# @examples
# data(bowhead.snps)
# snps <- df2gtypes(bowhead.snps, ploidy = 2)
# 
# snp.sfs <- sfs(snps, na.action = "resample")
# snp.sfs[1:10]
#' 
#' @name sfs
#' @export
#' 
sfs <- function(g, ref.allele = NULL, joint = FALSE, strata = NULL, 
                simplify = TRUE, na.action = c("fail", "filter", "resample")) {
  # check strata
  if(joint & !is.null(strata)) {
    missing <- setdiff(strata, getStrataNames(g))
    if(length(missing) > 0) {
      stop(
        "These strata cannot be found in 'g':", 
        paste(missing, collapse = ", ")
      )
    }
  }
  
  # convert g to numeric SNP matrix
  snp.mat <- as.data.frame(g, coded.snps = TRUE, ref.allele = ref.allele)
  
  if(!joint) {
    .sfs_1d(snp.mat, na.action) 
  } else {
    if(getNumStrata(g) == 1) {
      .sfs_1d(snp.mat, na.action) 
    } else {
      .sfs_2d(snp.mat, g, strata, simplify, na.action)
    }
  }
}


#' @noRd
#'
.refAlleleFreq <- function(snp.mat, na.action = c("fail", "filter", "resample")) {
  snp.mat <- if(any(is.na(snp.mat))) {
    # remove columns with all missing genotypes
    to.keep <- apply(snp.mat, 2, function(x) !all(is.na(x)))
    snp.mat <- snp.mat[, to.keep, drop = FALSE]
    na.action <- match.arg(na.action)
    if(na.action == "fail") stop("missing data present - cannot calculate SFS")
    if(na.action == "filter") { # remove columns with missing genotypes
      to.keep <- apply(snp.mat, 2, function(x) all(!is.na(x)))
      if(sum(to.keep) == 0) stop("all loci have missing data - cannot calculate SFS")
      snp.mat[, to.keep, drop = FALSE]
    } else { # select minimum sample size of non-missing genotypes from each column
      min.n <- min(apply(snp.mat, 2, function(x) sum(!is.na(x))))
      apply(snp.mat, 2, function(x) sample(x[!is.na(x)], min.n))
    }
  } else snp.mat
  colSums(snp.mat)
}


#' @noRd
#'
.sfs_1d <- function(snp.mat, na.action = c("fail", "filter", "resample")) {
  allele.freq <- .refAlleleFreq(snp.mat, na.action)
  # create data frame of site frequency spectrum
  sfs <- table(allele.freq = allele.freq) %>% 
    as.data.frame %>% 
    dplyr::mutate(
      allele.freq = as.numeric(as.character(.data$allele.freq))
    )
  # pad sfs with integer values with 0 count and return sfs
  zero.freqs <- setdiff(0:((nrow(snp.mat) * 2) + 1), sfs$allele.freq)
  if(length(zero.freqs) > 0) {
    sfs <- sfs %>% 
      dplyr::bind_rows(data.frame(allele.freq = zero.freqs, Freq = 0)) %>% 
      dplyr::arrange(.data$allele.freq)
  }
  stats::setNames(sfs$Freq, sfs$allele.freq)
}


#' @noRd
#'
.sfs_2d <- function(snp.mat, g, strata, simplify = TRUE,
                    na.action = c("fail", "filter", "resample")) {
  # split numerical SNP matrix into list of strata
  st.snp <- split(as.data.frame(snp.mat), getStrata(g)[rownames(snp.mat)])
  # use all strata if none specified
  if(is.null(strata)) strata <- getStrataNames(g)
  # get by-locus SNP frequencies for each stratum
  st.snp <- st.snp[strata]
  st.freqs <- do.call(cbind, lapply(st.snp, .refAlleleFreq, na.action = na.action))
  # get SFS for all pairs of strata
  st.pairs <- combn(strata, 2)
  result <- lapply(1:ncol(st.pairs), function(i) {
    st <- st.pairs[, i]
    # get count of all pairs of frequencies
    sfs <- outer(
      0:((nrow(st.snp[[st[1]]]) * 2) + 1),
      0:((nrow(st.snp[[st[2]]]) * 2) + 1),
      Vectorize(
        function(i, j, freqs) {
          sum(freqs[, 1] == i & freqs[, 2] == j)
        }, 
        c("i", "j")
      ),
      freqs = st.freqs[, st]
    )
    dimnames(sfs) <- stats::setNames(
      list(0:(nrow(sfs) - 1), 0:(ncol(sfs) - 1)),
      st
    )
    sfs
  })
  if(simplify & length(result) == 1) result[[1]] else result
}