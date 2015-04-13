####################
# Unions of classes (stolen from adegenet by Thibaut Jombart)
####################
setClassUnion("listOrNULL", c("list","NULL"))
setClassUnion("factorOrNULL", c("factor","NULL"))
setClassUnion("charOrNULL", c("character","NULL"))
setClassUnion("intOrNum", c("integer","numeric","NULL"))
setClassUnion("intOrNULL", c("integer","NULL"))
setClassUnion("dfOrNULL", c("data.frame", "NULL"))

#' @aliases gtypes
#' @import apex
#' @import methods
#'
#' @rdname gtypes
#'
#' @slot loci a data.frame containing the allelic data as one column per locus.
#'   Alleles are on multiple rows per column with samples listed in the same
#'   order for each allele. rownames are sample names plus allele number
#'   formatted as 12345.1 and 12345.2 where 12345 is the sample name and 1 and
#'   2 are the first and second alleles. colnames are unique locus names.
#' @slot sequences a \linkS4class{multidna} object, which is a list of
#'   \link{\code{DNAbin}} objects.
#' @slot schemes a data.frame with stratification schemes in each column.
#'   Sample names are in the rownames and must match the first part of the
#'   sample names (rownames) of the 'loci' slot. Each column is a factor.
#' @slot strata a factor as long as the number of samples representing the
#'   current stratification scheme.
#' @slot ploidy integer representing the ploidy of the data. There are
#'   ploidy * the number of samples rows in 'loci'.
#' @slot description a label for the object.
#' @slot other a slot to carry other related information - unused in package
#'   analyses.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#'
#' #--- create a microsatellite object
#' strata.schemes <- dolph.strata[, c("broad", "fine")]
#' rownames(strata.schemes) <- dolph.strata$id
#' msats <- gtypes.create(dolph.msats[, -1], ploidy = 2,
#'                        ind.names = dolph.msats[, 1],
#'                        schemes = strata.schemes)
#'
#'
#' #--- create a dloop object
#' dloop.haps <- dolph.strata$id
#' names(dloop.haps) <- dolph.strata$id
#' dloop.fine <- gtypes.create(dloop.haps, ploidy = 1, schemes = strata.schemes,
#'                             sequences = dolph.seqs, seq.names = "dLoop",
#'                             strata = "fine")
#'
#' @export

setClass(Class = "gtypes",
         slots = c(loci = "dfOrNULL", sequences = "dfOrNULL",
                   schemes = "dfOrNULL", strata = "factorOrNULL",
                   ploidy = "integer", description = "charOrNULL",
                   other = "ANY"
         ),
         prototype = prototype(loci = NULL, sequences = NULL, schemes = NULL,
                               strata = NULL, ploidy = 0, description = NULL,
                               other = NULL
         )
)


setValidity("gtypes", function(object) {
  # check that all columns in loci are factors
  loci.are.factors <- sapply(object@loci, is, class2 = "factor")
  if(!all(loci.are.factors)) {
    cat("all columns in the 'loci' slot are not factors\n")
    return(FALSE)
  }

  # check that length of sequences equals number of columns in loci
  num.seqs <- length(object@sequences@dna)
  if(num.seqs > 0 & num.seqs != ncol(object@loci)) {
    cat("the number of sets of sequences is not equal to the number of loci\n")
    return(FALSE)
  }

  # check that sequence haplotype labels can be found
  if(num.seqs > 0) {
    locus.good <- sapply(colnames(object@loci), function(x) {
      haps <- unique(as.character(object@loci[, x]))
      seqs <- rownames(object@sequences@dna[[x]])
      all(haps %in% seqs)
    })
    if(!all(locus.good)) {
      bad.loci <- paste(colnames(object@loci)[!locus.good], collapse = ", ")
      cat("haplotypes are missing in", bad.loci, "\n")
      return(FALSE)
    }
  }

  # check that ploidy is compatible with loci
  if(nrow(object@loci) %% object@ploidy != 0) {
    cat("number of alleles is not an even multiple of 'ploidy'\n")
    return(FALSE)
  }

  # check that strata is same length as number of individuals
  if(length(object@strata) > 0) {
    if((nrow(object@loci) / object@ploidy) != length(object@strata)) {
      cat("length of 'strata' slot is not equal to number of individuals\n")
      return(FALSE)
    }
  }

  # check that some individuals are in strata schemes
  if(nrow(object@schemes) > 0) {
    ids <- rownames(object@loci)[1:(nrow(object@loci) / object@ploidy)]
    ids <- substr(ids, 1, nchar(ids) - 2)
    if(length(intersect(ids, rownames(object@schemes))) == 0) {
      cat("no sample ids from the 'loci' slot are in stratification schemes\n")
      return(FALSE)
    }
  }

  TRUE
})
