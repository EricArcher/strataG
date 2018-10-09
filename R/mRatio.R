#' @title M ratio
#' @description Calculate Garza-Williamson M ratio (bottleneck) statistic for
#'   microsattelite data.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata calculate ratio for each stratum separately?
#' @param rpt.size set of values to check for allele repeat size. Function will
#'   use the largest common denominator found in this vector or return
#'   \code{NA}.
#'
#' @note The function will only compute the metric for microastellite loci,
#'   which is defined as loci with allele labels that can be converted to
#'   numeric values in their entirety and have a fixed repeat size. \code{NA} is
#'   returned for all loci that do not have all numeric alleles.  
#'   \code{NA} will also be returned if a locus is monomorphic, the locus has no
#'   genotypes, or a minimum repeat size cannot be found for all alleles at a
#'   locus.
#'   
#' @references Garza, J.C. and E.G. Williamson. 2001. Detection of reduction in
#'   population size using data form microsatellite loci. Molecular Ecology
#'   10(2):305-318.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' data(msats.g)
#' 
#' m.by.strata <- mRatio(msats.g, TRUE)
#' m.by.strata
#' 
#' m.overall <- mRatio(msats.g, FALSE)
#' m.overall
#' 
#' @export
#' 
mRatio <- function(g, by.strata = FALSE, rpt.size = 8:2) {
  # function takes a data.frame with allele and Freq
  .calcMratio <- function(df) {
    # check if all alleles are integers
    x <- sub("[[:space:]]+$", "", df$allele)
    x <- sub("^[[:space:]]+", "", x)
    x <- x[!x %in% c("", ".") | !is.na(x)]
    # no alleles are just character or numeric 
    if(!length(x)) return(NA)
    # some alleles can't be converted to integer
    if(suppressWarnings(any(is.na(as.integer(x))))) return(NA)
    # locus is monomorphic
    if(nrow(df) == 1) return(NA)
    # all frequencies are 0
    if(all(df$Freq == 0)) return(NA)
    
    # sort alleles in numerical order
    df <- df %>% 
      dplyr::mutate(allele = as.integer(.data$allele)) %>% 
      dplyr::arrange(.data$allele)
    
    # find repeat sizes
    size.diff <- diff(df$allele)
    rpt.found <- FALSE
    for(r in sort(rpt.size, decreasing = TRUE)) {
      if(all(size.diff %% r == 0)) {
        rpt.found <- TRUE
        break
      }
    }
    if(!rpt.found) return(NA)
    
    # number of alleles between smallest and largest that are present
    non.0 <- dplyr::filter(df, .data$Freq > 0) 
    n <- diff(non.0$allele[c(1, nrow(non.0))]) / r
    
    # calculate metric
    sum(df$Freq > 0) / (n + 1)
  }
  
  freqs <- alleleFreqs(g, by.strata = by.strata)
  freqs <- if(by.strata) {
    purrr::imap_dfr(freqs, function(f, x) {
      as.data.frame(f, stringsAsFactors = F) %>% 
        dplyr::rename(allele = .data$Var1, stratum = .data$Var2) %>% 
        dplyr::mutate(locus = x) 
    }) %>% 
      dplyr::group_by(.data$stratum, .data$locus)
  } else {
    purrr::imap_dfr(freqs, function(f, x) {
      as.data.frame(f, stringsAsFactors = F) %>% 
        dplyr::rename(allele = .data$Var1) %>% 
        dplyr::mutate(locus = x)
    }) %>% 
      dplyr::group_by(.data$locus)
  }

  freqs %>% 
    dplyr::do(m.ratio = .calcMratio(.data)) %>% 
    dplyr::ungroup() %>%
    tidyr::unnest(.data$m.ratio) %>% 
    as.data.frame()
}