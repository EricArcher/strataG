#' @title Input functions for fastsimcoal parameters
#' @description Functions to specify deme parameters, historical events, 
#'   migration matrices, and genetic data to simulate to write fastsimcoal 
#'   parameter or template files.
#'   
#' @param sequence.length number of base pairs to use for each block.
#' @param num.loci number of loci to simulate.
#' @param mut.rate per base pair or locus mutation rate.
#' @param transition.rate dna: fraction of substitutions that are transitions. 
#'   Set to 1 (all transitions) for SNPs.
#' @param gsm.param value of the geometric parameter for a
#'   Generalized Stepwise Mutation (GSM) model. This value represents the
#'   proportion of mutations that will change the allele size by more than
#'   one step. Values between 0 and 1 are required. A value of 0 is for a
#'   strict Stepwise Mutation Model (SMM).
#' @param range.constraint \code{msat}: Range constraint (number of different
#'   alleles allowed). A value of 0 means no range constraint.
#' @param recomb.rate recombination rate between adjacent markers. No effect for 
#'   SNPs.
#' @param chromosome number or character identifying which chromosome the marker
#'   is on.
#' @param deme.size the deme size.
#' @param sample.size the number of samples to take.
#' @param sample.time the number of generations in the past at which samples 
#'   are taken.
#' @param inbreeding the inbreeding coefficient for the deme \code{[0:1]}.
#' @param growth the growth rate of the deme.
#' @param event.time the number of generations before present at which the
#'   historical event happened.
#' @param source the source deme (the first listed deme has index 0).
#' @param sink the sink deme.
#' @param prop.migrants the expected proportion of migrants to move from 
#'   the source to the sink deme.
#' @param new.size the new size for the sink deme, relative to its size in 
#'   the previous (later in time) generation.
#' @param new.growth the new growth rate for the sink deme.
#' @param migr.mat the number of the new migration matrix to be used 
#'   further back in time. The matrices are those supplied to the 
#'   \code{fscSettingsMigration} function. The first matrix has index 0.
#' @param ploidy the desired ploidy of the final data. Ne and the number of 
#'   samples specified in \code{fscDeme} will be multiplied by this value as
#'   \code{fastsimcoal2} generates haploid data.
#' @param num.chrom the number of chromosomes to be simulated. If this is 
#'   specified and not the same as the number of linkage blocks specified by 
#'   the \code{fscMarker_} functions, then this many chromosomes with 
#'   duplicated structures will be simulated. If \code{NULL}, then
#'   the chromosome specification for each block will be used. 
#' @param output should the FREQ marker type be labelled for OUTPUT?
#' @param ... a set of comma-separated square migration matrices for 
#'   \code{fscSettingsMigration}, or marker specifications from a
#'   \code{fscMarker_} function for \code{fscSettingsGenetics}.
#'  
#' @note SNPs are simulated as a DNA sequence with a transiton rate of 1.
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
#' @seealso \code{\link{fscWrite}}
#' 
#' @name fsc.input
#' @aliases fscSettings
#' 
NULL


# ---- Demes ----

#' @rdname fsc.input
#' @export
#' 
fscDeme <- function(deme.size, sample.size, sample.time = 0, 
                    inbreeding = 0, growth = 0) {
  if(!(is.na(sample.size) | is.na(deme.size))) {
    if(sample.size > deme.size) {
      stop("`sample.size` can't be greater than `deme.size`.")
    }
  }
  
  deme <- c(
    deme.size = deme.size, sample.size = sample.size, sample.time = sample.time, 
    inbreeding = inbreeding, growth = growth
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
  demes <- do.call(rbind, demes)
  class(demes) <- c("fscSettingsDemes", class(demes))
  attr(demes, "ploidy") <- ploidy
  demes
}


# ---- Historical Events ----

#' @rdname fsc.input
#' @export
#' 
fscEvent <- function(event.time = 0, source = 0, sink = 0, prop.migrants = 1, 
                     new.size = 1, new.growth = 0, migr.mat = 0) {
  ev <- c(
    event.time = event.time, source = source, sink = sink, 
    prop.migrants = prop.migrants, new.size = new.size, 
    new.growth = new.growth, migr.mat = migr.mat
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
  if(all(!is.na(events[, 1]))) {
    events <- events[order(events[, 1]), , drop = FALSE]
  }
  class(events) <- c("fscSettingsEvents", class(events))
  events
}


# ---- Migration Rates ----

#' @rdname fsc.input
#' @export
#' 
fscSettingsMigration <- function(...) {
  migration <- list(...)
  for(mat in migration) {
    if(!(is.matrix(mat) & is.numeric(mat))) {
      stop("All values must be numeric matrices.")
    }
    if(nrow(mat) != ncol(mat)) stop("All matrices must be square.")
  }
  class(migration) <- c("fscSettingsMigration", class(migration))
  migration
}


# ---- Genetics ----

#' @rdname fsc.input
#' @export
#' 
fscMarker_dna <- function(sequence.length, mut.rate, recomb.rate = 0, 
                          transition.rate = 1/3, chromosome = 1) {
  loc <- data.frame(
    chrom = chromosome, 
    actual.type = "DNA",
    fsc.type = "DNA",
    num.markers = sequence.length,
    recomb.rate = recomb.rate,
    mut.rate = mut.rate,
    param.5 = transition.rate,
    param.6 = as.numeric(NA),
    stringsAsFactors = FALSE
  )
  class(loc) <- c("fscMarker", class(loc))
  loc
}

#' @name fsc.input
#' @export
#' 
fscMarker_microsat <- function(num.loci, mut.rate, recomb.rate = 0, 
                               gsm.param = 0, range.constraint = 0, 
                               chromosome = 1) {
  loc <- data.frame(
    chrom = chromosome,
    actual.type = "MICROSAT",
    fsc.type = "MICROSAT",
    num.markers = num.loci,
    recomb.rate = recomb.rate,
    mut.rate = mut.rate,
    param.5 = gsm.param,
    param.6 = range.constraint,
    stringsAsFactors = FALSE
  )
  class(loc) <- c("fscMarker", class(loc))
  loc
}

#' @rdname fsc.input
#' @export
#' 
fscMarker_snp <- function(sequence.length, mut.rate, 
                          recomb.rate = 0, chromosome = 1) {
  loc <- data.frame(
    chrom = chromosome, 
    actual.type = "SNP",
    fsc.type = "DNA",
    num.markers = sequence.length,
    recomb.rate = recomb.rate,
    mut.rate = mut.rate,
    param.5 = 1,
    param.6 = as.numeric(NA),
    stringsAsFactors = FALSE
  )
  class(loc) <- c("fscMarker", class(loc))
  loc
}

#' @rdname fsc.input
#' @export
#' 
fscMarker_standard <- function(num.loci, mut.rate, recomb.rate = 0,
                               chromosome = 1) {
  loc <- data.frame(
    chrom = chromosome,
    actual.type = "STANDARD",
    fsc.type = "STANDARD",
    num.markers = num.loci,
    recomb.rate = recomb.rate,
    mut.rate = mut.rate,
    param.5 = as.numeric(NA),
    param.6 = as.numeric(NA),
    stringsAsFactors = FALSE
  )
  class(loc) <- c("fscMarker", class(loc))
  loc
}

#' @rdname fsc.input
#' @export
#'
fscMarker_freq <- function(mut.rate, output = TRUE) {  
  loc <- data.frame(
    chrom = 1,
    actual.type = "FREQ",
    fsc.type = "FREQ",
    num.markers = 1,
    recomb.rate = 0,
    mut.rate = mut.rate,
    param.5 = ifelse(output, "OUTEXP", ""),
    param.6 = as.numeric(NA),
    stringsAsFactors = FALSE
  )
  class(loc) <- c("fscMarker", class(loc))
  loc
}

#' @rdname fsc.input
#' @export
#' 
fscSettingsGenetics <- function(..., num.chrom = NULL) {
  if(!is.null(num.chrom)) {
    if(!is.numeric(num.chrom)) stop("`num.chrom` must be numeric.")
    if(num.chrom < 1) stop("`num.chrom` can't be < 1.")
  }
  
  markers <- list(...)
  
  # check supplied markers
  for(x in markers) {
    if(!inherits(x, "fscMarker")) {
      stop("Marker definitions must be produced by a `fscMarker_xxx()` function.")
    }
    if(x$fsc.type == "FREQ" & length(markers) > 1) {
      stop("Marker type FREQ can't be specified with other markers.")
    }
  }
  
  # identify chromosome structure
  chrom <- do.call(rbind, markers)
  if(is.null(num.chrom)) {
    chrom <- lapply(split(chrom, chrom$chrom), function(x) x[, -1])
    num.chrom <- length(chrom)
  } else chrom <- list(chrom[, -1])
  chrom.same <- num.chrom == 1 | (num.chrom > 1 & length(chrom) == 1)
  
  attr(chrom, "num.chrom") <- num.chrom
  attr(chrom, "chrom.diff") <- !chrom.same
  class(chrom) <- c("fscSettingsGenetics", class(chrom))
  chrom
}