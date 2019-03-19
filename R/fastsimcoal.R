#' @title Run fastsimcoal
#' @description Run a fastsimcoal simulation and load results into a 
#'   \linkS4class{gtypes} object.
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
#' @param infinite.alleles use infinite alleles model?
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
#' @seealso \code{\link{fastsimcoal.input}}
#' 
#' @name fastsimcoal
#' 
NULL


# ---- Writing ----

#' @noRd
.fscWritePar <- function(p) {
  fname <- paste0(p$label, ifelse(p$is.tpl, ".tpl", ".par"))
  f <- file(fname, open = "wt")
  
  n.demes <- nrow(p$demes)
  writeLines("//Number of population samples (demes)", f)
  writeLines(as.character(n.demes), f)
  
  writeLines("//Population effective sizes (number of genes)", f)
  for(i in 1:n.demes) writeLines(p$demes[i, 1], f)
  
  writeLines("//Sample sizes, times, inbreeding", f)
  for(i in 1:n.demes) writeLines(paste(p$demes[i, 2:4], collapse = " "), f)
  
  writeLines("//Growth rates: negative growth implies population expansion", f)
  for(i in 1:n.demes) writeLines(p$demes[i, 5], f)
  
  writeLines("//Number of migration matrices: 0 implies no migration between demes", f)
  writeLines(as.character(length(p$migration)), f)
  if(!is.null(p$migration)) {
    for(i in 1:length(p$migration)) {
      writeLines("//migration matrix", f)
      for(r in 1:nrow(p$migration[[i]])) {
        writeLines(paste(p$migration[[i]][r, ], collapse = " "), f)
      }
    }
  }
  
  writeLines("//Historical events: time, source, sink, migrants, new size, growth rate, migr. matrix", f)
  n.events <- if(is.null(p$events)) 0 else nrow(p$events)
  writeLines(as.character(n.events), f)
  if(!is.null(p$events)) {
    for(r in 1:n.events) writeLines(paste(p$events[r, ], collapse = " "), f)
  }
  
  writeLines("//Number of independent loci [chromosomes]", f)
  num.chrom <- attributes(p$genetics)["num.chrom"]
  chrom.diff <- as.numeric(attributes(p$genetics)["chrom.diff"])
  writeLines(paste(num.chrom, chrom.diff, collapse = " "), f)
  for(block in p$genetics) {
    writeLines("//Per chromosome: Number of linkage blocks", f)
    writeLines(as.character(nrow(block)), f)
    writeLines("//Per block: data type, num loci, rec. rate and mut rate + optional parameters", f)
    for(i in 1:nrow(block)) {
      writeLines(paste(as.character(block[i, -1]), collapse = " "), f)
    }
  }
  
  close(f)
  p$in.file <- fname
  p
}


#' @noRd
.fscWriteEst <- function(p) {
  fname <- paste0(p$label, ".est")
  p$est.file <- fname
  
  if(file.exists(fname) & interactive()) {
    prompt <- paste0(
      "The estimation file, '", 
      fname, 
      "' already exists. Overwrite? (Y/n) "
    )
    if(readline(prompt) != "Y") return(p)
  }
  
  f <- file(fname, open = "wt")
  
  est.df <- data.frame()
  for(col in colnames(p$demes)) {
    param <- grep("[$]", p$demes[, col], value = T)
    if(length(param) > 0) {
      is.int <- col %in% c("deme.size", "sample.size", "sample.time")
      est.df <- rbind(est.df, data.frame(int = as.numeric(is.int), param))
    }
  }
  
  if(!is.null(p$events)) {
    for(col in colnames(p$events)) {
      param <- grep("[$]", p$events[, col], value = T)
      if(length(param) > 0) {
        is.int <- col %in% c("event.time", "source", "sink", "migr.mat")
        est.df <- rbind(est.df, data.frame(int = as.numeric(is.int), param))
      }
    }
  }
  
  if(!is.null(p$migration)) {
    for(mat in p$migration) {
      param <- grep("[$]", mat, value = T)
      if(length(param) > 0) est.df <- rbind(est.df, data.frame(int = 0, param))
    }
  }  
  
  for(chrom in p$genetics) {
    for(col in colnames(chrom)) {
      param <- grep("[$]", chrom[, col], value = T)
      if(length(param) > 0) {
        is.int <- col == "param.6"
        est.df <- rbind(est.df, data.frame(int = as.numeric(is.int), param))
      }
    }
  }
  
  colnames(est.df) <- c("int", "param")
  est.df$param <- as.character(est.df$param)
  
  writeLines("// Priors and rules file", f)
  writeLines("// *********************", f)
  writeLines("", f)
  writeLines("[PARAMETERS]", f)
  writeLines("//#isInt? #name #dist. #min #max", f)
  writeLines("//all N are in number of haploid individuals", f)
  for(i in 1:nrow(est.df)) {
    writeLines(
      paste0(
        paste(est.df[i, ], collapse = "\t"),
        "\tunif\tmin\tmax\toutput", 
        collapse = ""
      ), 
      f
    )
  }
  writeLines("", f)
  writeLines("[RULES]", f)
  writeLines("", f)
  writeLines("[COMPLEX PARAMETERS]", f)
  for(i in 1:nrow(est.df)) {
    writeLines(
      paste0(paste(est.df[i, ], collapse = "\t"), " = ", collapse = ""), 
      f
    )
  }
  
  close(f)
  
  message(
    "The estimation file, '", 
    fname, 
    "', has been written and must be edited prior to running."
  )
  
  invisible(p)
}


#' @rdname fastsimcoal
#' @export
#' 
fscWrite <- function(demes, genetics, migration = NULL, events = NULL, 
                     label = "strataG.fastsimcoal") {
  opt <- options(scipen = 999)
  params <- list(
    demes = demes,
    migration = migration, 
    events = events,
    genetics = genetics
  )
  params <- .checkFscInput(params) 
  params$label <- make.names(label)
  params <- .fscWritePar(params) 
  if(params$is.tpl) params <- .fscWriteEst(params) 
  options(opt)
  
  invisible(params)
}


# ---- Running ----

#' @rdname fastsimcoal
#' @export
#' 
fscRun <- function(p, num.sims = 1, dna.to.snp = FALSE, max.snps = NULL, 
                   all.sites = TRUE, infinite.alleles = FALSE, 
                   sfs.type = c("daf", "maf"), num.ecm.loops = 20,
                   num.cores = NULL, seed = NULL, exec = "fsc26") {
  
  if(file.exists(p$label)) unlink(p$label, recursive = TRUE, force = TRUE) 
  
  opt <- options(scipen = 999)
  
  cores.spec <- if(!is.null(num.cores)) {
    num.cores <- max(1, num.cores)
    num.cores <- min(num.cores, min(parallel::detectCores(), 12))
    if(num.cores == 1) "" else {
      paste(c("--cores", "--numBatches"), num.cores, collapse = " ")
    }
  } else ""
  
  if(!is.null(max.snps)) dna.to.snp <- TRUE
  sfs.type <- switch(match.arg(sfs.type), daf = "--dsfs", maf = "--msfs")
  
  args <- c(
    ifelse(p$is.tpl, "--tplfile", "--ifile"), p$in.file,
    "--numsims", num.sims
  )
  if(p$is.tpl) args <- c(args, "-e", p$est.file)
  if(!is.null(seed)) args <- c(args, "--seed ", seed)
  if(dna.to.snp) args <- c(
    args, "--dnatosnp", ifelse(is.null(max.snps), 0, max.snps)
  )
  if(all.sites) args <- c(args, "--allsites")
  if(infinite.alleles) args <- c(args, "--inf")
  if(p$is.tpl) args <- c(args, sfs.type)
  if(p$is.tpl) args <- c(args, "--maxlhood")
  if(p$is.tpl) args <- c(args, "--numloops", num.ecm.loops)
  args <- c(args, cores.spec)
  args <- paste(args, collapse = " ")
  
  p$run.params <- list(
    num.sims = num.sims, dna.to.snp = dna.to.snp, max.snps = max.snps, 
    all.sites = all.sites, infinite.alleles = infinite.alleles, 
    sfs.type = sfs.type, num.ecm.loops = num.ecm.loops,
    num.cores = num.cores, seed = seed, exec = exec
  )
  p$args <- args
  p$log.file <- paste0(p$label, ".log")
  if(file.exists(p$log.file)) file.remove(p$log.file)
  cat(format(Sys.time()), "running fastsimcoal...\n")
  err <- system2(exec, p$args, stdout = p$log.file)
  # err <- if(.Platform$OS.type == "unix") {
  #   system2(exec, args, stdout = p$log.file)
  # } else {
  #   shell(paste(exec, args), intern = F)
  # }
  if(err != 0) {
    stop(
      format(Sys.time()), 
      "fastsimcoal exited with error ", err, "\n",
      "The command was:\n",
      exec, p$args
    )
  }
  
  Sys.sleep(1) # wait for .arp files to finish writing
  p$arp.files <- NULL
  arp.files <- dir(p$label, pattern = ".arp$", full.names = TRUE)
  if(length(arp.files) > 0) {
    arp.files <- arp.files[order(nchar(arp.files), arp.files)]
    p$arp.files <- stats::setNames(
      arp.files, 
      paste0("rep", 1:length(arp.files))
    )
    cat(format(Sys.time()), "creating locus map for .arp files...\n")
    p <- .fscMapArpLocusInfo(p)
  }
  
  options(opt)
  cat(format(Sys.time()), "run complete\n")
  invisible(p)
}


# ---- Reading ----

#' @noRd
.zeroPad <- function(x) {
  formatC(x, digits = floor(log10(max(x))), flag = "0", mode = "integer") 
}

#' @noRd
.fscMapArpLocusInfo <- function(p) {
  # expand genetic info in input parameters to matrix
  locus.info <- do.call(rbind, p$genetics)
  if(!attr(p$genetics, "chrom.diff")) {
    num.blocks <- nrow(locus.info)
    num.chrom <- attr(p$genetics, "num.chrom")
    locus.info <- replicate(num.chrom, locus.info, simplify = FALSE)
    locus.info <- do.call(rbind, locus.info)
    locus.info <- cbind(
      chromosome = rep(1:num.chrom, each = num.blocks),
      locus.info
    )
  }
  
  # create list to associate rows in locus info to columns read from .arp file
  .extractMatCol <- function(i, next.col, locus.info) {
    if(locus.info$fsc.type[i] == "DNA") { # DNA takes 1 column
      next.col
    } else { # MICROSAT and STANDARD take num.markers columns
      next.col:(next.col + locus.info$num.markers[i] - 1)
    }
  }
  mat.col <- vector("list", nrow(locus.info)) 
  mat.col[[1]] <- .extractMatCol(1, 3, locus.info)
  if(nrow(locus.info) > 1) {
    for(i in 2:nrow(locus.info)) {
      # repeat column number for consecutive DNA entries
      if(locus.info$fsc.type[i] == "DNA" & locus.info$fsc.type[i - 1] == "DNA") {
        mat.col[[i]] <- mat.col[[i - 1]]
      } else {
        next.col <- max(mat.col[[i - 1]]) + 1
        mat.col[[i]] <- .extractMatCol(i, next.col, locus.info)
      }
    }
  }
  
  # add column of locus identifiers to locus information
  locus.info <- locus.info %>% 
    dplyr::mutate(mat.col = mat.col) %>% 
    dplyr::group_by(.data$chromosome) %>% 
    dplyr::mutate(block = 1:dplyr::n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      name = paste0(
        "C", .zeroPad(.data$chromosome),
        "_B", .zeroPad(.data$block), 
        "_", .data$actual.type
      ),
      dna.col = NA
    ) %>% 
    dplyr::ungroup()
  
  # identify columns of consecutive DNA blocks
  for(i in 1:nrow(locus.info)) {
    if(locus.info$fsc.type[i] == "DNA") {
      locus.info[i, "dna.col"] <- locus.info$mat.col[[i]]
    }
  }
  
  # identify start and stop characters for each DNA block within its column
  p$locus.info <- locus.info %>% 
    dplyr::group_by(.data$dna.col) %>% 
    dplyr::mutate(
      dna.end = ifelse(
        .data$fsc.type == "DNA", 
        cumsum(.data$num.markers), 
        NA
      ),
      dna.start = ifelse(
        .data$fsc.type == "DNA", 
        .data$dna.end - .data$num.markers + 1, 
        NA
      )
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(.data$name, .data$chromosome, .data$block, dplyr::everything()) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  
  invisible(p)
}


#' @rdname fastsimcoal
#' @export
#' 
fscReadArpFile <- function(file) {
  # read .arp file
  f <- readr::read_lines(file)
  
  # get start and end points of data blocks
  start <- grep("SampleData=", f) + 1
  end <- which(f == "}") - 2
  pos <- cbind(start, end)
  
  # extract matrix for each data block
  data.mat <- do.call(rbind, lapply(1:nrow(pos), function(i, pos) {
    f.line <- f[pos[i, 1]:pos[i, 2]]
    result <- do.call(rbind, strsplit(f.line, "[[:space:]]+"))[, -2]
    cbind(result[, 1], paste0("Deme.", rep(i, nrow(result))), result[, -1])
  }, pos = pos))
  
  colnames(data.mat) <- c("id", "deme", paste0("col", 3:ncol(data.mat)))
  data.mat
}


#' @noRd
#' 
.fscParseLocusRows <- compiler::cmpfun(function(locus.info, data.mat, .zeroPad) {
  # for each row in locus.info, create matrix of loci in that block
  # return list of block matrices
  locus.rows <- vector("list", nrow(locus.info))
  for(i in 1:nrow(locus.info)) {
    cols <- data.mat[, locus.info$mat.col[[i]]]
    # extract DNA sequence from string in column
    if(locus.info$fsc.type[i] == "DNA") {
      cols <- stringi::stri_sub(
        cols, 
        locus.info$dna.start[i], 
        locus.info$dna.end[i]
      )
      if(locus.info$actual.type[i] == "SNP") {
        cols <- do.call(rbind, strsplit(cols, ""))
      }
    }
    cols <- cbind(cols)
    # give suffixes to block names with more than one locus
    colnames(cols) <- if(ncol(cols) == 1) {
      locus.info$name[i] 
    } else {
      paste0(locus.info$name[i], "_L", .zeroPad(1:ncol(cols)))
    }
    locus.rows[[i]] <- cols
  }
  locus.rows
})


#' @rdname fastsimcoal
#' @export
#' 
fscParseGeneticData <- function(p, sim = 1) {
  if(any(p$locus.info$fsc.type == "DNA") & !p$run.params$all.sites) {
    stop(
      "Can't read .arp output if DNA sequences are present and fastsimcoal has not been run with `all.sites = TRUE`\n"
    )
  }
  if(!is.numeric(sim) | length(sim) > 1) stop("`sim` must be a single number")
  
  file <- p$arp.files[sim]
  cat(format(Sys.time()), "reading", file, "\n")
  data.mat <- fscReadArpFile(file)
  
  # parse data matrix based on entries in locus.info data.frame
  cat(format(Sys.time()), "parsing locus data...\n")
  gen.data <- .fscParseLocusRows(p$locus.info, data.mat, .zeroPad)
  
  # create map of column numbers in gen.data for each row of locus.info
  last.col <- 2
  locus.cols <- vector("list", nrow(p$locus.info))
  names(locus.cols) <- p$locus.info$name
  for(i in 1:length(gen.data)) {
    locus.cols[[i]] <- (last.col + 1):(last.col + ncol(gen.data[[i]]))
    last.col <- max(locus.cols[[i]])
  }
  
  gen.data <- do.call(cbind, gen.data)
  gen.data <- cbind(data.mat[, 1:2], gen.data)  
  attr(gen.data, "locus.cols") <- locus.cols
  attr(gen.data, "file") <- file
  
  invisible(gen.data)
}


#' @rdname fastsimcoal
#' @export
#' 
fscReadVector <- function(file) {    
  f <- scan(file, what = "character", sep = "\n", quiet = TRUE)
  stats::setNames(
    as.numeric(do.call(c, strsplit(f[2], "\t"))),
    do.call(c, strsplit(f[1], "\t"))
  )
}


#' @rdname fastsimcoal
#' @export
#' 
fscReadEstParams <- function(p) {
  sfs.file <- dir(
    p$label, 
    pattern = paste0("^", p$label, "_[[:alnum:]]+.txt$"),
    full.names = TRUE
  )  
  brent.file <- dir(
    p$label, 
    pattern = paste0("^", p$label, ".brent_lhoods$"),
    full.names = TRUE
  )
  best.file <- dir(
    p$label, 
    pattern = paste0("^", p$label, ".bestlhoods$"),
    full.names = TRUE
  )
  if(length(sfs.file) == 0 | length(brent.file) == 0 | length(best.file) == 0) {
    stop("Can't file all output files (*.txt, *.brent_lhoods, *.bestlhoods)")
  }
  
  brent.file <- readr::read_lines(brent.file) 
  brent.file <- brent.file[grep("^Param|^[[:digit:]]", brent.file)]
  brent.file <- strsplit(brent.file, "\t")
  brent <- as.data.frame(do.call(rbind, lapply(brent.file[-1], as.numeric)))
  colnames(brent) <- brent.file[[1]][1:ncol(brent)]
  
  sfs.vec <- if(length(sfs.file) == 1) {
    fscReadVector(sfs.file)
  } else {
    sapply(sfs.file, fscReadVector, simplify = FALSE)
  }
  
  invisible(
    list(
      ml.params = fscReadVector(best.file),
      brent.lhoods = brent,
      sfs = sfs.vec
    )
  )
}


#' @rdname fastsimcoal
#' @export
#' 
fscRead <- function(p, sim = 1) {
  if(p$is.tpl) fscReadEstParams(p) else fscParseGeneticData(p, sim)
}



#' @rdname fastsimcoal
#' @export
#' 
fscExtractLoci <- function(p, sim = 1, gen.data = NULL, type = "all", 
                           sep.chrom = FALSE, chrom = NULL) {
  if(is.null(gen.data)) {
    if(p$is.tpl) {
      stop("Can't extract loci because 'p' specifies a parameter estimation model.")
    }
    gen.data <- fscRead(p, sim)
  } 
  file <- attr(gen.data, "file")
  
  # filter locus info for specified chromosomes
  locus.info <- p$locus.info
  if(!is.null(chrom)) {
    if(!is.numeric(chrom)) stop("'chrom' must be a numeric vector")
    if(max(chrom) > max(locus.info$chromosome)) {
      stop("there are not", max(chrom), "chromosomes available") 
    }
    locus.info <- locus.info[locus.info$chromosome %in% chrom, ]
  }
  
  # check marker type
  type <- toupper(type)
  if(!all(type %in% c("DNA", "SNP", "MICROSAT", "STANDARD", "ALL"))) {
    stop("`type` can only contain 'dna', 'snp', 'microsat', 'standard', 'all'.")
  }
  if("ALL" %in% type) type <- unique(locus.info$actual.type)
  
  # filter locus info for specified marker types
  locus.info <- locus.info[grep(paste(type, collapse = "|"), locus.info$name), ]
  
  .extractLocCols <- function(loc.names, mat) {
    loc.cols <- unlist(attr(mat, "locus.cols")[loc.names]) 
    mat[, c(1:2, loc.cols), drop = FALSE]
  }
  
  # extract columns for selected chromosomes and marker types
  gen.data <- if(sep.chrom) { # extract for each chromosome
    locus.info$chrom.label <- regmatches(
      locus.info$name,
      regexpr("^C[[:alnum:]]+", locus.info$name)
    )
    tapply( 
      locus.info$name, 
      locus.info$chrom.label,
      .extractLocCols,
      mat = gen.data
    )
  } else { # extract selected loci from data frame
    .extractLocCols(locus.info$name, gen.data)
  }
  
  attr(gen.data, "file") <- file
  gen.data
}


#' @rdname fastsimcoal
#' @export
#' 
fscCleanup <- function(label, save.est = FALSE) {
  # remove label folder
  unlink(label, recursive = TRUE, force = TRUE)
  
  files <- c(dir(pattern = paste0("^", label)), "seed.txt")
  if(length(files) != 0) {
    # don't remove R script files
    r.files <- grep("[[:alnum:]]+\\.r$", files, ignore.case = TRUE, value = TRUE)
    if(length(r.files) != 0) files <- setdiff(files, r.files)
    if(save.est) { # leave .est files if they are to be saved
      est.files <- grep("[[:alnum:]]+\\.est$", files, ignore.case = TRUE, value = TRUE)
      if(length(est.files) != 0) files <- setdiff(files, est.files)
    }
    # remove only files, not other directories that start with label
    files <- files[utils::file_test("-f", files)]
    file.remove(files)
    invisible(files)
  } else invisible(NULL)
}


# ---- gtypes ----


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


#' @rdname fastsimcoal
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
