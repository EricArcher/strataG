#' @title Load fastsimcoal output to gtypes
#' @description Create a gtypes object from the result of a fastsimcoal run.
#'
#' @param demes matrix of deme sampling information created by the 
#'   \code{\link{fscSettingsDemes}} function.
#' @param genetics data.frame specifying loci to simulate created by the 
#'   \code{\link{fscSettingsGenetics}} function.
#' @param migration a list of matrices giving the migration rates 
#'   between pairs of demes created by the \code{\link{fscSettingsMigration}} 
#'   function.
#' @param events matrix of historical events created by the 
#'   \code{\link{fscSettingsEvents}} function.
#' @param label character string to label files with.
#' @param seed random number seed for simulation.
#' @param exec name of fastsimcoal executable.
#' @param num.cores number of cores to use.
#' @param file filename to read from or write to.
#' @param p list of fastsimcoal input parameters and output.
#' @param num.sims number of simulation replicates to run.
#' @param dna.to.snp convert DNA sequences to numerical SNPs?
#' @param max.snps maximum number of SNPs to retain.
#' @param all.sites retain all sites?
#' @param infinite.sites use infinite alleles model?
#' @param sfs.type type of site frequency spectrum to compute for each 
#'   population sample: `daf` = derived allele frequency (unfolded), 
#'   `maf` = minor allele frequency (folded).
#' @param num.ecm.loops number of loops (ECM cycles) to be performed when 
#'   estimating parameters from SFS. Default is 20.
#' @param save.est do not delete .est parameter estimation files during cleanup?
#' @param sim number of the simulation replicate to read.
#' @param gen.data matrix of parsed genetic data read from .arp file with 
#'   `fscParseGeneticData()`.
#' @param type type of marker to return.
#' @param sep.chrom return a list with chromosomes separated?
#' @param chrom numerical vector giving chromosomes to return. `NULL` 
#'   returns all chromosomes.
#' @param drop.mono drop monomorphic sites before creating `gtypes` object?
#' 
#' @note fastsimcoal is not included with `strataG` and must be downloaded 
#'   separately. Additionally, it must be installed such that it can be run from 
#'   the command line in the current working directory. See the vignette 
#'   for \code{external.programs} for installation instructions.
#' 
#' @references Excoffier, L. and Foll, M (2011) fastsimcoal: a continuous-time 
#'   coalescent simulator of genomic diversity under arbitrarily complex 
#'   evolutionary scenarios Bioinformatics 27: 1332-1334.\cr
#'   Excoffier, L., Dupanloup, I., Huerta-Sánchez, E., Sousa, V.C., 
#'   and M. Foll (2013) Robust demographic inference from genomic and SNP data. 
#'   PLOS Genetics, 9(10):e1003905. \cr
#'   \url{http://cmpg.unibe.ch/software/fastsimcoal2/}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{fsc.input}}
#' 
#' @export
#' 
fsc2gtypes <- function(p, sim = 1, gen.data = NULL, type = NULL, chrom = NULL, 
                       drop.mono = TRUE) {
  # check arguments
  p.types <- tolower(unique(p$locus.info$actual.type))
  type <- if(is.null(type)) {
    if(length(p.types) == 1) p.types else {
      stop(
        "fastsimcoal output must have a single locus type.",
        paste("Select one of the following types:", paste(p.types, collapse = ", "))
      )
    }
  } else {
    if(length(type) > 1 | !is.character(type)) {
      stop("`type` must be a single character string.")
    } else if(!type %in% p.types) {
      stop(
        "`type` was not one of the available locus types:",
        paste(p.types, collapse = ", ")
      )
    }
  }
  
  ploidy <- attr(p$genetics, "ploidy")
  if(ploidy == 1 & type != "dna") {
    stop("Can't create a haploid `gtypes` object for `type` other than 'dna'.")
  }
  sep.chrom <- ifelse(type == "dna" & ploidy == 1, TRUE, FALSE)
  
  # extract requested data if necessary
  if(is.null(gen.data)) {
    gen.data <- fscExtractLoci(p, sim, gen.data, type, sep.chrom, chrom)
  }
  description <- attr(gen.data, "file")
  
  # remove monomorphic loci
  if(drop.mono & ploidy > 1) {
    not.mono <- apply(gen.data[, -(1:2), drop = FALSE], 2, function(loc) {
      dplyr::n_distinct(loc) > 1
    })
    gen.data <- gen.data[, c(1:2, which(not.mono) + 2), drop = FALSE]
    if(ncol(gen.data) <= 2) {
      warning("All loci are monomorphic. NULL returned.")
      return(NULL)
    }
  }
  
  # create gtypes object
  cat(format(Sys.time()), "loading gtypes object...\n")
  if(ploidy > 1) { # compile non-haploid genotypes
    loc.names <- colnames(gen.data)[-(1:2)]
    loc.names <- paste(rep(loc.names, each = ploidy), 1:ploidy, sep = ".")
    gen.data <- t(sapply(seq(1, nrow(gen.data), by = ploidy), function(i) {
      rows <- i:(i + ploidy - 1)
      id <- paste(gen.data[rows, "id"], collapse = "|")
      c(id, gen.data[i, "deme"], as.vector(gen.data[rows, -(1:2)]))
    }))
    colnames(gen.data) <- c("id", "deme", loc.names)
    as.data.frame(gen.data, stringsAsFactors = FALSE) %>% 
      .loadGtypes(ploidy = ploidy, description = description)
  } else { # create list of sequences
    dna.list <- sapply(gen.data, function(chrom.mat) {
      apply(
        as.matrix(chrom.mat[, -(1:2)]),
        1,
        paste,
        collapse = ""
      ) %>% 
        strsplit(split = "") %>% 
        stats::setNames(chrom.mat[, "id"]) %>% 
        ape::as.DNAbin()
    }, simplify = FALSE, USE.NAMES = TRUE)
    gen.data <- gen.data[[1]][, c(1, 2, rep(1, length(dna.list))), drop = FALSE]
    colnames(gen.data)[3:ncol(gen.data)] <- names(dna.list)
    .loadGtypes(gen.data, ploidy, dna.list, description) %>% 
      labelHaplotypes()
  }
}

#' @noRd
#' 
.loadGtypes <- function(gen.mat, ploidy, sequences = NULL, description = NULL) {
  # return new gtypes object
  methods::new(
    "gtypes", 
    gen.data = gen.mat[, -(1:2)], 
    ploidy = ploidy, 
    ind.names = gen.mat[, 1],
    strata = gen.mat[, 2], 
    schemes = NULL, 
    sequences = sequences, 
    description = description, 
    other = list()
  )
}