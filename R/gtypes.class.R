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
#'   \code{\link{genind2gtypes}}, \code{\link{gtypes.accessors}},
#'   \code{\link{initialize.gtypes}}
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
#' dloop.g <- labelHaplotypes(dloop.g, "Hap.")$gtypes
#' dloop.g
#' 
#' @aliases gtypes
#' @import adegenet ape apex data.table
#' @importFrom methods setClass
#' @export
#' 
setClass(
  Class = "gtypes",
  slots = c(
    data = "data.table", sequences = "dnaSequences",
    ploidy = "integer", schemes = "data.frameOrNULL", 
    description = "charOrNULL", other = "ANY"
  ),
  prototype = prototype(
    data = NULL, sequences = NULL, ploidy = 0L, 
    schemes = NULL, description = NULL, other = NULL
  ),
  validity = function(object) {
    # check that the first two columns are "ids" and "strata"
    if(!identical(colnames(object@data)[1:2], c("ids", "strata"))) {
      cat("first two columns in the 'data' slot must be 'ids' and 'strata'\n")
      return(FALSE)
    }
    
    # check that all columns in data are factors
    loci.are.factors <- object@data[, sapply(.SD, is.factor), .SDcols = !c("ids", "strata")]
    if(!all(loci.are.factors)) {
      cat("all locus columns in the 'data' slot are not factors\n")
      return(FALSE)
    }

    # check sequences
    if(!is.null(object@sequences)) {
      # check that length of sequences equals number of locus columns
      dna <- getSequences(sequences(object), simplify = FALSE)
      num.seqs <- length(dna)
      if(num.seqs > 0 & num.seqs != ncol(object@data) - 2) {
        cat("the number of sets of sequences is not equal to the number of loci\n")
        return(FALSE)
      }

      # check that locus names are the same in the @data colnames and names of
      #  @sequences
      loc.names <- colnames(object@data)[-(1:2)]
      if(!identical(loc.names, getLocusNames(object@sequences))) {
        cat("the names of the sets of sequences is not the same as the loci\n")
        return(FALSE)
      }

      # check that sequence haplotype labels can be found
      locus.good <- sapply(loc.names, function(x) {
        haps <- object@data[[x]]
        seqs <- rownames(as.matrix(dna[[x]]))
        all(na.omit(haps) %in% seqs)
      })
      if(!all(locus.good)) {
        bad.loci <- paste(loc.names[!locus.good], collapse = ", ")
        cat("haplotypes are missing in:", bad.loci, "\n")
        return(FALSE)
      }
    }

    # check that ploidy is compatible with loci
    if(nrow(object@data) %% object@ploidy != 0) {
      cat("number of alleles is not an even multiple of 'ploidy'\n")
      return(FALSE)
    }
    if(object@ploidy != 1 & !is.null(object@sequences)) {
      cat("sequences can't be present unless object is haploid (ploidy = 1)\n")
      return(FALSE)
    }

    # check that at least some individuals are in strata schemes
    if(!is.null(object@schemes)) {
      if(length(intersect(object@data$ids, rownames(object@schemes))) == 0) {
        cat("no sample ids from the 'data' slot are in stratification schemes\n")
        return(FALSE)
      }
    }

    TRUE
  }
)