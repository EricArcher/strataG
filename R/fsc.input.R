#' @title Input functions for fastsimcoal parameters
#' @description These functions specify and format simulation parameters used to
#'   write fastsimcoal2 parameter or template files, parameter estimation files, 
#'   parameter definition files, and site frequency spectrum files.
#'   
#' @param sequence.length number of base pairs to use for each block.
#' @param num.loci number of loci to simulate.
#' @param mut.rate per base pair or locus mutation rate.
#' @param transition.rate dna: fraction of substitutions that are transitions.
#' @param gsm.param value of the geometric parameter for a Generalized Stepwise
#'   Mutation (GSM) model. This value represents the proportion of mutations
#'   that will change the allele size by more than one step. Values between 0
#'   and 1 are required. A value of 0 is for a strict Stepwise Mutation Model
#'   (SMM).
#' @param range.constraint \code{msat}: Range constraint (number of different
#'   alleles allowed). A value of 0 means no range constraint.
#' @param recomb.rate recombination rate between adjacent markers. No effect for
#'   SNPs.
#' @param chromosome number or character identifying which chromosome the marker
#'   is on.
#' @param deme.size the number of individuals in the deme.
#' @param sample.size the number of samples to take.
#' @param sample.time the number of generations in the past at which samples are
#'   taken.
#' @param inbreeding the inbreeding coefficient for the deme \code{[0:1]}.
#' @param growth the growth rate of the deme.
#' @param event.time the number of generations before present at which the
#'   historical event happened.
#' @param source the source deme (the first listed deme has index 0).
#' @param sink the sink deme.
#' @param prop.migrants the expected proportion of migrants to move from the
#'   source to the sink deme.
#' @param new.size the new size for the sink deme, relative to its size in the
#'   previous (later in time) generation.
#' @param new.growth the new growth rate for the sink deme.
#' @param migr.mat the number of the new migration matrix to be used further
#'   back in time. The matrices are those supplied to the
#'   \code{fscSettingsMigration} function. The first matrix has index 0.
#' @param ploidy the desired ploidy of the final data. \code{deme.size} and 
#'   \code{sample.size} will be multiplied by this value in the parameter or 
#'   template file as \code{fastsimcoal2} generates haploid data.
#' @param num.chrom the number of chromosomes to be simulated. If this is
#'   specified and not the same as the number of linkage blocks specified by the
#'   \code{fscBlock_} functions, then this many chromosomes with duplicated
#'   structures will be simulated. If \code{num.chrom = NULL}, then the
#'   chromosome specification for each block will be used.
#' @param outexp logical describing if the expected site frequency
#'   spectrum given the estimated parameters should be output?
#' @param name name of the parameter being specified. Must match a name used in
#'   one of the simulation settings functions.
#' @param is.int logical specifying whether or not the parameter is an integer.
#' @param distr a character string giving the distribution to use to select initial values for
#'   parameter estimation. Can be \code{"unif"} or \code{"logunif"}.
#' @param min,max minimum and maximum values for the distribution specified in \code{distr}.
#' @param value character string giving the value that the complex parameter is
#'   to take.
#' @param output logical indicating if estimates for the parameter should be
#'   output.
#' @param bounded logical indicating whether to treat the parameter as a bounded
#'   estimate.
#' @param reference logical indicating whether the parameter is to be used as a
#'   reference.
#' @param obs.sfs vector, matrix, or list containing observed SFS to use for
#'   parameter estimation.
#' @param rules character vector giving rules for the parameter estimation.
#' @param sfs.type type of SFS to write. Can be \code{maf} or \code{daf}.
#' @param mat numeric matrix or data frame with values of parameters to use in
#'   place of parameter names in simulation.
#' @param ... a set of comma-separated values for settings. See 
#'   Notes for more information.
#'
#' @note All settings must be passed to \code{\link{fscWrite}} using one
#'   of the \code{fscSettingsXXX} functions. Most of these functions in turn
#'   take as their input comma-separated values which are the result of specific
#'   \code{fscXXX} functions:\cr
#'   \describe{ 
#'     \item{fscSettingsDemes()}{comma-separated instances of fscDeme(). If 
#'        names are given for each deme, these names will be used in the parsed 
#'        output.} 
#'     \item{fscSettingsEvents()}{comma-separated instances of fscEvent().} 
#'     \item{fscSettingsMigration()}{comma-separated migration matrices.} 
#'     \item{fscSettingsGenetics()}{comma-separated instances of 
#'       fscBlock_dna(), fscBlock_microsat(), fscBlock_snp(), 
#'       fscBlock_standard(), or fscBlock_freq(). SNPs are simulated as a 
#'       DNA sequence with a transiton rate of 1. `fscBlock_freq()` can only 
#'       be used by itself and in parameter estimation simulations.} 
#'     \item{fscSettingsEst()}{comma-separated instances of fscEstParam() 
#'       as well as site frequency spectra (\code{obs.sfs}) and
#'       parameter rules {\code{rules}}.}
#'   }
#'   \code{fastsimcoal2} is not included with `strataG` and must be downloaded 
#'     separately. Additionally, it must be installed such that it can be run
#'     from the command line in the current working directory. 
#'     The function \code{fscTutorial()} will open a detailed tutorial on the 
#'     interface in your web browser.
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
#' @seealso \code{\link{fscWrite}}, \code{\link{fscRun}}, \code{\link{fscRead}}
#'  
#' @examples
#' # three demes with optional names
#' demes <- fscSettingsDemes(
#'   Large = fscDeme(10000, 10), 
#'   Small = fscDeme(2500, 10),
#'   Medium = fscDeme(5000, 3, 1500)
#' )
#' 
#' # four historic events
#' events <- fscSettingsEvents(
#'   fscEvent(event.time = 2000, source = 1, sink = 2, prop.migrants = 0.05),
#'   fscEvent(2980, 1, 1, 0, 0.04),
#'   fscEvent(3000, 1, 0),
#'   fscEvent(15000, 0, 2, new.size = 3)
#'  )
#'  
#' # four genetic blocks of different types on three chromosomes.  
#' genetics <- fscSettingsGenetics(
#'   fscBlock_snp(10, 1e-6, chromosome = 1),
#'   fscBlock_dna(10, 1e-5, chromosome = 1),
#'   fscBlock_microsat(3, 1e-4, chromosome = 2),
#'   fscBlock_standard(5, 1e-3, chromosome = 3)
#' )
#' 
#' #' same four genetic blocks of different types with same structure repeated three times.  
#' genetics <- fscSettingsGenetics(
#'   fscBlock_snp(10, 1e-6),
#'   fscBlock_dna(10, 1e-5),
#'   fscBlock_microsat(3, 1e-4),
#'   fscBlock_standard(5, 1e-3),
#'   num.chrom = 3
#' )
#'
#' @name fsc.input
#' @aliases fscSettings
#' 
NULL

#' @noRd
#' 
.checkSingleArgs <- function() {
  for(x in ls(envir = parent.frame())) {
    x.val <- eval(parse(text = x), envir = parent.frame())
    if(!is.vector(x.val)) stop("'", x, "' must be a vector.", call. = FALSE)
    if(length(x.val) > 1) stop("'", x, "' must be one value.", call. = FALSE)
  }
}


# Demes ------------------------------------------------------------------------

#' @rdname fsc.input
#' @export
#' 
fscDeme <- function(deme.size, sample.size, sample.time = 0, 
                    inbreeding = 0, growth = 0) {
  .checkSingleArgs()
  deme <- data.frame(
    deme.size = deme.size, sample.size = sample.size, 
    sample.time = sample.time, inbreeding = inbreeding, growth = growth,
    stringsAsFactors = FALSE
  )
  class(deme) <- c("fscDeme", class(deme))
  deme
}


#' @rdname fsc.input
#' @export
#' 
fscSettingsDemes <- function(..., ploidy = 2) {
  demes <- list(...)
  for(x in demes) if(!inherits(x, "fscDeme")) {
    stop("All demes must be produced by `fscDeme()`.")
  }
  if(is.null(names(demes))) {
    names(demes) <- paste0("Deme.", 1:length(demes))
  } else {
    i <- which(names(demes) == "")
    names(demes)[i] <- paste0("Deme.", i)
  }
  demes <- do.call(rbind, demes) %>% 
    tibble::rownames_to_column("deme.name")
  class(demes) <- c("fscSettingsDemes", class(demes))
  attr(demes, "ploidy") <- ploidy
  demes
}


# Historical Events ------------------------------------------------------------

#' @rdname fsc.input
#' @export
#' 
fscEvent <- function(event.time = 0, source = 0, sink = 0, prop.migrants = 1, 
                     new.size = 1, new.growth = 0, migr.mat = 0) {
  .checkSingleArgs()
  ev <- data.frame(
    event.time = event.time, source = source, sink = sink, 
    prop.migrants = prop.migrants, new.size = new.size, 
    new.growth = new.growth, migr.mat = migr.mat,
    stringsAsFactors = FALSE
  )
  class(ev) <- c("fscEvent", class(ev))
  ev
}

#' @rdname fsc.input
#' @export
#' 
fscSettingsEvents <- function(...) {
  events <- list(...)
  for(x in events) if(!inherits(x, "fscEvent")) {
    stop("Historical events must be produced by `fscEvent()`.")
  }
  events <- do.call(rbind, events)
  class(events) <- c("fscSettingsEvents", class(events))
  events
}


# Migration Rates --------------------------------------------------------------

#' @rdname fsc.input
#' @export
#' 
fscSettingsMigration <- function(...) {
  migration <- list(...)
  for(mat in migration) {
    if(!is.matrix(mat)) stop("All values must be matrices.")
    if(nrow(mat) != ncol(mat)) stop("All matrices must be square.")
  }
  class(migration) <- c("fscSettingsMigration", class(migration))
  migration
}


# Genetics ---------------------------------------------------------------------

.makeBlock <- function(chromosome, actual.type, fsc.type, num.markers,
                       recomb.rate, mut.rate, param.5, param.6) {
  block <- data.frame(
    chromosome = as.integer(chromosome), 
    actual.type = actual.type,
    fsc.type = fsc.type,
    num.markers = as.integer(num.markers),
    recomb.rate = recomb.rate,
    mut.rate = mut.rate,
    param.5 = param.5,
    param.6 = param.6,
    stringsAsFactors = FALSE
  )
  class(block) <- c("fscBlock", class(block))
  block
}

#' @rdname fsc.input
#' @export
#' 
fscBlock_dna <- function(sequence.length, mut.rate, recomb.rate = 0, 
                         transition.rate = 1/3, chromosome = 1) {
  .checkSingleArgs()
  if(!is.numeric(sequence.length)) stop("'sequence.length' must be a number.")
  .makeBlock(
    chromosome = chromosome, actual.type = "DNA", fsc.type = "DNA",
    num.markers = sequence.length, recomb.rate = recomb.rate,
    mut.rate = mut.rate, param.5 = transition.rate, param.6 = as.numeric(NA)
  )
}

#' @name fsc.input
#' @export
#' 
fscBlock_microsat <- function(num.loci, mut.rate, recomb.rate = 0, 
                              gsm.param = 0, range.constraint = 0, 
                              chromosome = 1) {
  .checkSingleArgs()  
  if(!is.numeric(num.loci)) stop("'num.loci' must be a number.")
  .makeBlock(
    chromosome = chromosome, actual.type = "MICROSAT", fsc.type = "MICROSAT",
    num.markers = num.loci, recomb.rate = recomb.rate, mut.rate = mut.rate,
    param.5 = gsm.param, param.6 = range.constraint
  )
}

#' @rdname fsc.input
#' @export
#' 
fscBlock_snp <- function(sequence.length, mut.rate, recomb.rate = 0, 
                         chromosome = 1) {
  .checkSingleArgs()
  if(!is.numeric(sequence.length)) stop("'sequence.length' must be a number.")
  .makeBlock(
    chromosome = chromosome, actual.type = "SNP", fsc.type = "DNA",
    num.markers = sequence.length, recomb.rate = recomb.rate,
    mut.rate = mut.rate, param.5 = 1, param.6 = as.numeric(NA)
  )
}

#' @rdname fsc.input
#' @export
#' 
fscBlock_standard <- function(num.loci, mut.rate, recomb.rate = 0,
                              chromosome = 1) {
  .checkSingleArgs()
  if(!is.numeric(num.loci)) stop("'num.loci' must be a number.")
  .makeBlock(
    chromosome = chromosome, actual.type = "STANDARD", fsc.type = "STANDARD",
    num.markers = num.loci, recomb.rate = recomb.rate, mut.rate = mut.rate,
    param.5 = as.numeric(NA), param.6 = as.numeric(NA)
  )
}

#' @rdname fsc.input
#' @export
#'
fscBlock_freq <- function(mut.rate, outexp = TRUE) {  
  .checkSingleArgs()
  .makeBlock(
    chromosome = 1L, actual.type = "FREQ", fsc.type = "FREQ", num.markers = 1L,
    recomb.rate = 0L, mut.rate = mut.rate,
    param.5 = ifelse(outexp, "OUTEXP", ""), param.6 = as.numeric(NA)
  )
}

#' @rdname fsc.input
#' @export
#' 
fscSettingsGenetics <- function(..., num.chrom = NULL) {
  if(!is.null(num.chrom)) {
    if(!is.numeric(num.chrom)) stop("`num.chrom` must be numeric.")
    if(num.chrom < 1) stop("`num.chrom` can't be < 1.")
  }
  
  blocks <- list(...)
  
  # check supplied markers
  for(x in blocks) {
    if(!inherits(x, "fscBlock")) {
      stop("Block definitions must be produced by `fscBlock_xxx()` functions.")
    }
    if(x$fsc.type == "FREQ" & length(blocks) > 1) {
      stop("Block type FREQ can't be specified with other blocks")
    }
  }
  
  blocks <- do.call(rbind, blocks)
  
  # identify chromosome structure
  chrom.same <- TRUE
  if(is.null(num.chrom)) {
    num.chrom <- length(unique(blocks$chromosome))
    if(num.chrom > 1) chrom.same <- FALSE
  }
  
  attr(blocks, "num.chrom") <- num.chrom
  attr(blocks, "chrom.diff") <- !chrom.same
  class(blocks) <- c("fscSettingsGenetics", class(blocks))
  blocks
}


# Parameter Estimation ---------------------------------------------------------

#' @rdname fsc.input
#' @export
#' 
fscEstParam <- function(name, is.int = TRUE, distr = c("unif", "logunif"), 
                        min = NA, max = NA, value = NA, output = TRUE, 
                        bounded = FALSE, reference = FALSE) {
  distr <- match.arg(distr)
  .checkSingleArgs()
  if(is.na(value) & (is.na(min) | is.na(max))) {
    stop("Either 'min' and 'max' or 'value' must be specified.")
  }
  param.df <- data.frame(
    is.int = as.integer(is.int),
    name = toupper(name),
    dist = ifelse(is.na(value), distr, as.character(NA)),
    min = ifelse(is.na(value), as.character(min), as.character(NA)),
    max = ifelse(is.na(value), as.character(max), as.character(NA)),
    value = ifelse(!is.na(value), as.character(value), as.character(NA)),
    output = ifelse(output, "output", "hide"),
    bounded = ifelse(bounded, "bounded", ""),
    reference = ifelse(reference, "reference", ""),
    stringsAsFactors = FALSE
  )
  class(param.df) <- c("fscEstParam", class(param.df))
  param.df
}

#' @rdname fsc.input
#' @export
#' 
fscSettingsEst <- function(..., obs.sfs, rules = NULL, 
                           sfs.type = c("maf", "daf")) {
  params <- list(...)
  for(x in params) {
    if(!"fscEstParam" %in% class(x)) {
      stop("All parameters in '...' must be generated by 'fscEstParam()'.")
    }
  }
  params <- do.call(rbind, params)
  
  if(!(is.vector(obs.sfs) | is.matrix(obs.sfs) | is.list(obs.sfs))) {
    stop("'obs.sfs' must be a vector, matrix, or list")
  }
  
  if(!is.null(rules) & !is.character(rules)) {
    stop("'rules' must be a character vector.")
  }
  
  if(is.vector(obs.sfs) & !is.list(obs.sfs)) {
    if(!is.numeric(obs.sfs)) stop("'obs.sfs' must be a numeric vector.")
    names(obs.sfs) <- paste0("d0_", 0:(length(obs.sfs) - 1))
  }
  
  if(is.matrix(obs.sfs)) {
    if(!is.numeric(obs.sfs)) stop("'obs.sfs' must be a numeric matrix.")
    rownames(obs.sfs) <- paste0("d1_", 0:(nrow(obs.sfs) - 1))
    colnames(obs.sfs) <- paste0("d0_", 0:(ncol(obs.sfs) - 1))
    obs.sfs <- list(obs.sfs)
  }
  
  if(is.list(obs.sfs)) {
    .fscFormattedNames <- function(x.names) {
      if(is.null(x.names)) return(FALSE)
      all(grepl("^d[[:digit:]]+_[[:digit:]]+$", x.names))
    }
    
    for(mat in obs.sfs) {
      good.rownames <- .fscFormattedNames(rownames(mat))
      good.colnames <- .fscFormattedNames(colnames(mat))
      if(!good.rownames | !good.colnames) {
        stop(
          "All matrices in 'obs.sfs' must have fastsimcoal formatted ",
          "rowmames and colnames (e.g., 'd0_0', 'd0_1', 'd0_2', etc)."
        )
      }
    }
  }
  
  est <- list(params = params, rules = rules, sfs = obs.sfs)
  attr(est, "sfs.type") <- toupper(match.arg(sfs.type))
  class(est) <- c("fscSettingsEst", class(est))
  est
}


# Parameter Definition ---------------------------------------------------------

#' @rdname fsc.input
#' @export
#' 
fscSettingsDef <- function(mat) {
  if(is.data.frame(mat)) mat <- as.matrix(mat)
  if(!is.numeric(mat) & !is.matrix(mat)) stop("'mat' must be a numeric matrix.")
  if(is.null(colnames(mat))) stop("'mat' must have colnames.")
  class(mat) <- c("fscSettingsDef", class(mat))
  mat
}
