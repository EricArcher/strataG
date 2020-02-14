setClassUnion("dnaSequences", c("multidna", "NULL"))

#' @title \code{gtypes} Class
#' @description An S4 class storing multi-allelic locus or sequence data along
#'   with a current stratification and option stratification schemes.
#'
#' @slot data a data.table where the first column contains the sample ID 
#'   (\code{ids}). The second column contains the sample stratification 
#'   (\code{strata}). The third column to the end contains the allelic data as 
#'   one column per locus. Alleles are on multiple rows per column with sample 
#'   IDs duplicated for all alleles. Column names are unique locus names.
#' @slot sequences a \linkS4class{multidna} object.
#' @slot ploidy integer representing the ploidy of the data. There are
#'   ploidy * the number of individuals rows in 'data'.
#' @slot schemes a data.frame with stratification schemes in each column.
#'   The rownames are individual names and must match the 'id' column of the
#'   'data' slot. Each column is a factor.
#' @slot description a label for the object (optional).
#' @slot other a slot to carry other related information - currently unused in 
#'   analyses (optional).
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @seealso \code{\link{df2gtypes}}, \code{\link{sequence2gtypes}},
#'   \code{\link{gtypes.accessors}}, \code{\link{gtypes.initialize}}
#' @examples
#'
#' #--- create a diploid (microsatellite) gtypes object
#' data(dolph.msats)
#' data(dolph.strata)
#' strata.schemes <- dolph.strata[, c("broad", "fine")]
#' rownames(strata.schemes) <- dolph.strata$id
#' msats.g <- new("gtypes", gen.data = dolph.msats[, -1], ploidy = 2,
#'                ind.names = dolph.msats[, 1], schemes = strata.schemes)
#' msats.g
#'
#' #--- create a haploid sequence (mtDNA) gtypes object and label haplotypes
#' data(dolph.seqs)
#' dloop.haps <- cbind(dLoop = dolph.strata$id)
#' rownames(dloop.haps) <- dolph.strata$id
#' dloop.g <- new("gtypes", gen.data = dloop.haps, ploidy = 1, 
#'                schemes = strata.schemes, sequences = dolph.seqs, 
#'                strata = "fine")
#' dloop.g
#' labelHaplotypes(dloop.g, "Hap.")
#' 
#' @aliases gtypes
#' @import data.table apex
#' @export
#' 
methods::setClass(
  Class = "gtypes",
  slots = c(
    data = "data.table", sequences = "dnaSequences",
    ploidy = "integer", schemes = "data.frameOrNULL", 
    description = "charOrNULL", other = "list"
  ),
  
  prototype = prototype(
    data = NULL, 
    sequences = NULL, 
    ploidy = 0L, 
    schemes = NULL, 
    description = NULL, 
    other = list()
  ),
  
  validity = function(object) {
    data.cols <- colnames(object@data)
    
    # check that there are only four columns in 'data'
    if(length(data.cols) != 4) {
      cat("'data' slot does not have 4 columns\n")
      return(FALSE)
    }
    
    # check that columns "id", "stratum", "locus", and "allele" are present
    if(!all(data.cols %in% c("id", "stratum", "locus", "allele"))) {
      cat("column names of 'data' slot must be 'id', 'stratum', 'locus', and 'allele'\n")
      return(FALSE)
    }
    
    # check that all columns in data are character
    if(!all(sapply(object@data, is.character))) {
      cat("all columns in the 'data' slot are not characters\n")
      return(FALSE)
    }

    # check sequences
    if(!is.null(object@sequences)) {
      # check that length of sequences equals number of loci
      dna <- getSequences(object)
      num.seqs <- length(dna)
      if(num.seqs > 0 & num.seqs != length(unique(object@data[["locus"]]))) {
        cat("the number of sets of sequences is not equal to the number of loci\n")
        return(FALSE)
      }

      # check that locus names are the same in the @data 'locus' column and
      #   names of @sequences
      loc.names <- sort(unique(object@data[["locus"]]))
      if(!identical(loc.names, sort(apex::getLocusNames(object@sequences)))) {
        cat("the names of the sets of sequences is not the same as the loci\n")
        return(FALSE)
      }

      # check that sequence haplotype labels can be found
      locus.good <- sapply(loc.names, function(x) {
        haps <- object@data %>% 
          dplyr::filter(.data$locus == x) %>% 
          dplyr::pull(.data$allele)
        seqs <- rownames(as.matrix(dna[[x]]))
        all(na.omit(haps) %in% seqs)
      })
      if(!all(locus.good)) {
        bad.loci <- paste(loc.names[!locus.good], collapse = ", ")
        cat("haplotypes are missing in:", bad.loci, "\n")
        return(FALSE)
      }
    }

    # check that ploidy is the same for all individuals/loci
    pl <- unique(table(object@data[["id"]], object@data[["locus"]]))
    if(length(pl) != 1) {
      cat("some individuals have different numbers of alleles per locus\n")
      return(FALSE)
    }
    
    # check that true ploidy matches stored ploidy
    if(object@ploidy != pl) {
      cat("ploidy in 'data' slot not same as in 'ploidy' slot\n")
      return(FALSE)
    }
    
    # check that sequences aren't present if not haploid
    if(object@ploidy != 1 & !is.null(object@sequences)) {
      cat("sequences can't be present unless object is haploid (ploidy = 1)\n")
      return(FALSE)
    }

    # check that at least some individuals are in strata schemes
    if(!is.null(object@schemes)) {
      if(length(intersect(object@data$id, object@schemes$id)) == 0) {
        cat("no sample ids from the 'data' slot are in stratification schemes\n")
        return(FALSE)
      }
    }

    TRUE
  }
)