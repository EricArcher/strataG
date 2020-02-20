#' @name as.data.frame.gtypes
#' @title Convert \code{gtypes} to data.frame or matrix
#' @description Create a formatted data.frame or matrix 
#'   from a \linkS4class{gtypes} object.
#'   
#' @param x a \linkS4class{gtypes} object.
#' @param one.col logical. If \code{TRUE}, then result has one column per 
#'   locus.
#' @param sep character to use to separate alleles if \code{one.col} is 
#'   \code{TRUE}.
#' @param ids logical. include a column for individual identifiers (\code{ids})?
#' @param strata logical. include a column for current statification 
#'   (\code{strata})?
#' @param sort.alleles logical. for non-haploid objects, should alleles be
#'   sorted in genotypes or left in original order? (only takes affect if
#'   \code{one.col = TRUE})
#' @param coded.snps return diploid SNPs coded as 0 (reference allele homozygote), 
#'   1 (heterozygote), or 2 (alternate allele homozygote). If this is `TRUE`,  
#'   the data is diploid, and all loci are biallelic, a data frame of 
#'   coded genotypes will be returned with one column per locus.
#' @param ref.allele an optional vector of reference alleles for each SNP. 
#'   Only used if `coded.snps = TRUE`. If provided, it must be at least as long 
#'   as there are biallelic SNPs in \code{g}. If named, the names must 
#'   match those of all loci in \code{g}. If set to `NULL` (default) the 
#'   major allele at each SNP is used as the reference.
#' @param ... additional arguments to be passed to or from methods.
#'   
#' @return A \code{data.frame} or \code{matrix} with one row per individual.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{df2gtypes} \link[strataG]{as.matrix}
#' 
#' @examples 
#' data(msats.g)
#' 
#' # with defaults (alleles in multiple columns, with ids and stratification)
#' df <- as.data.frame(msats.g)
#' str(df)
#' 
#' # one column per locus
#' onecol.df <- as.data.frame(msats.g, one.col = TRUE)
#' str(onecol.df)
#' 
#' # just the genotypes
#' genotypes.df <- as.data.frame(msats.g, one.col = TRUE, ids = FALSE, strata = FALSE)
#' str(genotypes.df)
#' 
#' # as a matrix instead
#' genotypes.mat <- as.matrix(msats.g)
#' str(genotypes.mat)
#' 
#' @aliases as.data.frame,gtypes-method as.data.frame.gtypes as.data.frame
#'
#' @export
#' 
methods::setMethod(
  "as.data.frame", 
  "gtypes",
  function(x, one.col = FALSE, sep = "/", ids = TRUE, 
           strata = TRUE, sort.alleles = TRUE, coded.snps = FALSE,
           ref.allele = NULL, ...) {
    
    if(coded.snps) { 
      if(getPloidy(x) != 2) stop("Can't code SNPs in non-diploid data.")
      # check that all loci are biallelic
      if(!all(numAlleles(x)$num.alleles <= 2)) {
        stop("Can't code SNPs because some loci have more than 2 alleles.")
      }
      
      # check reference allele vector
      if(!is.null(ref.allele)) {
        if(is.null(names(ref.allele))) {
          if(length(ref.allele) != getNumLoci(x)) {
            stop(
              "`ref.allele` is not as long as the number of loci and",
              "is not named."
            )
          } 
          names(ref.allele) <- getLociNames(x)
        } else {
          missing.loci <- setdiff(getLociNames(x), names(ref.allele))
          if(length(missing.loci) > 0) {
            stop(
              "The following loci can't be found in `ref.allele`:",
              paste(missing.loci, collapse = ", ")
            )
          }
        }
      } else {
        # set reference to major allele if ref.allele not supplied
        af <- alleleFreqs(x)
        ref.allele <- sapply(alleleFreqs(x), function(loc) {
          names(which.max(loc))
        })
      }
      
      # matrix of genotypes with one column per locus
      x.df <- x@data %>% 
        dplyr::mutate(ref = ref.allele[.data$locus]) %>% 
        dplyr::group_by(.data$id, .data$stratum, .data$locus) %>% 
        dplyr::summarize(code = sum(.data$allele == .data$ref)) %>% 
        dplyr::ungroup() %>% 
        tidyr::spread(.data$locus, .data$code) %>% 
        as.data.frame
      
      # remove ids or strata if requested
      if(!ids) x.df$id <- NULL
      if(!strata) x.df$stratum <- NULL
      
      x.df
    } else {
      # create data.frame of one column per locus
      x.df <- x@data %>% 
        dplyr::group_by(.data$id, .data$stratum, .data$locus) %>% 
        dplyr::summarize(
          genotype = .combineLoci(.data$allele, sep = sep, sort = sort.alleles)
        ) %>% 
        dplyr::ungroup() %>% 
        tidyr::spread("locus", "genotype") 
      
      if(getPloidy(x) == 1) one.col <- TRUE
      # if loci are to be split into separate columns, use alleleSplit
      if(!one.col) {
        x.df <- cbind(
          x.df[, c("id", "stratum")],
          x.df %>% 
            dplyr::select(-.data$id, -.data$stratum) %>% 
            as.data.frame() %>% 
            alleleSplit(sep = sep),
          stringsAsFactors = FALSE
        )
      }
      
      x.df <- dplyr::select(x.df, .data$id, .data$stratum, dplyr::everything())
      
      # remove ids or strata if requested
      if(!ids) x.df$id <- NULL
      if(!strata) x.df$stratum <- NULL
      
      x.df %>% 
        dplyr::mutate_all(dplyr::funs("as.character")) %>% 
        as.data.frame()
    }
  })


#' @rdname as.data.frame.gtypes
#' @aliases as.matrix,gtypes-method as.matrix.gtypes as.matrix
#' 
#' @export
#' 
methods::setMethod(
  "as.matrix", 
  "gtypes",
  function(x, one.col = FALSE, sep = "/", ids = TRUE, strata = TRUE, 
           sort.alleles = TRUE,...) {
    as.matrix(as.data.frame(
      x = x,
      one.col = one.col,
      sep = sep,
      ids = ids,
      strata = strata,
      sort.alleles = sort.alleles,
      ...
    ))
  })
