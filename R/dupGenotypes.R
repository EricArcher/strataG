#' @title Duplicate Genotypes
#' @description Identify duplicate or very similar genotypes.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param num.shared either number of loci or percentage of loci two 
#'   individuals must share to be considered duplicate individuals.
#' 
#' @return if no duplicates are present, the result is \code{NULL}, otherwise
#'   a data frame with the following columns is returned:
#' \tabular{ll}{
#'   \code{ids.1, ids.2} \tab sample ids.\cr
#'   \code{strata.1, strata.2} \tab sample stratification.\cr
#'   \code{mismatch.loci} \tab loci where the two samples do not match.\cr
#'   \code{num.loci.genotyped} \tab number of loci genotyped for both 
#'     samples.\cr
#'   \code{num.loci.shared} \tab number of loci shared (all alleles the same) between both samples.\cr
#'   \code{prop.loci.shared} \tab proportion of loci genotyped for both samples 
#'     that are shared.\cr
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' 
#' # identify potential duplicates in Coastal strata
#' coastal <- msats.g[, , "Coastal"]
#' coastal.5 <- coastal[getIndNames(coastal)[1:5], , ]
#' dupes <- dupGenotypes(coastal.5)
#' dupes
#' 
#' @export
#' 
dupGenotypes <- function(g, num.shared = 0.8) {
  #if not already, convert num.shared to %
  if(num.shared > 1) num.shared <- num.shared / getNumLoci(g) 
    
  dup.df <- propSharedLoci(g, type = "ids") %>% 
    dplyr::filter(.data$prop.loci.shared >= num.shared)
  
  dup.df <- if(nrow(dup.df) > 0) {
    st <- getStrata(g)
    locs <- getLociNames(g)
    
    dup.df %>% 
      tidyr::gather(
        "locus", "prop.shared", 
        -(.data$ids.1:.data$prop.loci.shared)
      ) %>% 
      dplyr::filter(.data$prop.shared < 1) %>% 
      dplyr::group_by(.data$ids.1, .data$ids.2) %>% 
      dplyr::summarize(mismatch.loci = paste(.data$locus, collapse = ", ")) %>% 
      dplyr::ungroup() %>% 
      dplyr::right_join(
        dplyr::select(dup.df, .data$ids.1:.data$prop.loci.shared),
        by = c("ids.1", "ids.2")
      ) %>% 
      dplyr::mutate(
        strata.1 = as.character(st[.data$ids.1]),
        strata.2 = as.character(st[.data$ids.2])
      ) %>% 
      dplyr::arrange(
        dplyr::desc(.data$prop.loci.shared), dplyr::desc(.data$num.loci.shared),
        dplyr::desc(.data$ids.1), dplyr::desc(.data$ids.2)
      ) %>% 
      dplyr::select(
        .data$ids.1, .data$ids.2, 
        .data$strata.1, .data$strata.2, dplyr::everything()
      ) %>% 
      as.data.frame()
  } else NULL
  
  if(is.null(dup.df)) cat("No duplicates found. NULL returned.\n")
  dup.df
}
