#' @title Duplicate Genotypes
#' @description Identify duplicate or very similar genotypes.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param num.shared either number of loci or percentage of loci two 
#'   individuals must share to be considered duplicate individuals.
#' @param num.cores number of CPU cores to use.
#' 
#' @return if no duplicates are present, the result is \code{NULL}, otherwise
#'   a data.frame with the following columns is returned:
#' \tabular{ll}{
#'   \code{ids.1, ids.2} \tab sample ids.\cr
#'   \code{strata.1, strata.2} \tab sample stratification.\cr
#'   \code{num.loci.genotyped} \tab number of loci genotyped for both 
#'     samples.\cr
#'   \code{num.loci.shared} \tab number of loci shared (all alleles the same) between both samples.\cr
#'   \code{prop.loci.shared} \tab proportion of loci genotyped for both samples 
#'     that are shared.\cr
#'   \code{mismatch.loci} \tab loci where the two samples do not match.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' 
#' # identify potential duplicates in Coastal strata
#' dupes <- dupGenotypes(msats.g[, , "Coastal"])
#' dupes
#' 
#' @export
#' 
dupGenotypes <- function(g, num.shared = 0.8, num.cores = NULL) {
  #if not already, convert num.shared to %
  if(num.shared > 1) num.shared <- num.shared / getNumLoci(g) 
    
  dup.df <- propSharedLoci(g, type = "ids", num.cores = num.cores) %>% 
    dplyr::filter(prop.loci.shared >= num.shared)
  
  dup.df <- if(nrow(dup.df) > 0) {
    st <- strata(g)
    locs <- getLocusNames(g)
    dup.df %>% 
      tidyr::gather(locus, prop.shared, -(ids.1:prop.loci.shared)) %>% 
      dplyr::filter(prop.shared < 1) %>% 
      dplyr::group_by(ids.1, ids.2) %>% 
      dplyr::summarize(mismatch.loci = paste(locus, collapse = ", ")) %>% 
      dplyr::ungroup() %>% 
      dplyr::right_join(
        dplyr::select(dup.df, ids.1:prop.loci.shared),
        by = c("ids.1", "ids.2")
      ) %>% 
      dplyr::mutate(
        strata.1 = as.character(st[ids.1]),
        strata.2 = as.character(st[ids.2])
      ) %>% 
      dplyr::arrange(
        dplyr::desc(prop.loci.shared),
        dplyr::desc(num.loci.shared),
        dplyr::desc(ids.1),
        dplyr::desc(ids.2)
      ) %>% 
      dplyr::select(ids.1, ids.2, strata.1, strata.2, dplyr::everything()) %>% 
      as.data.frame()
  } else NULL
  
  if(is.null(dup.df)) cat("No duplicates found. NULL returned.\n")
  dup.df
}
