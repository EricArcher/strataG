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
dupGenotypes <- function(g, num.shared = 0.8, num.cores = 1) {
  #if not already, convert num.shared to %
  if(num.shared > 1) num.shared <- num.shared / nLoc(g) 
    
  shared.locs <- propSharedLoci(g, type = "ids", num.cores = num.cores)
  dup.df <- shared.locs[shared.locs[, "prop.same"] >= num.shared, ]
  if(nrow(dup.df) > 0) {
    dup.df$strata.1 <- as.character(strata(g)[dup.df$ids.1])
    dup.df$strata.2 <- as.character(strata(g)[dup.df$ids.2])
    dup.df$mismatch.loci <- sapply(1:nrow(dup.df), function(i) {
      shared.prop <- as.matrix(dup.df[i, locNames(g)])
      loc.diff <- locNames(g)[which(shared.prop < 1)]
      paste(loc.diff, collapse = ", ")
    })
    colnames(dup.df)[c(3:5)] <- c(
      "num.loci.shared", "num.loci.genotyped", "prop.loci.shared"
    )
    dup.df <- dup.df[, c("ids.1", "ids.2", "strata.1", "strata.2", 
                         "num.loci.genotyped", "num.loci.shared", 
                         "prop.loci.shared", "mismatch.loci")]
  } 
  
  if(nrow(dup.df) > 0) {
    sort.order <- order(dup.df$prop.loci.shared, dup.df$num.loci.shared, 
                        rev(dup.df$ids.1), rev(dup.df$ids.2), decreasing = TRUE
    )
    dup.df <- dup.df[sort.order, ]
    rownames(dup.df) <- NULL
  } else dup.df <- NULL
  
  if(is.null(dup.df)) cat("No duplicates found. NULL returned.\n")
  dup.df
}
