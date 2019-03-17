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
#' @param delete.files logical. Delete files when done?
#' @param read.output logical. Read output of simulation when done?
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
                   num.cores = NULL, seed = NULL, delete.files = FALSE, 
                   read.output = TRUE, exec = "fsc26") {
  
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
  
  p$log.file <- paste0(p$label, ".log")
  if(file.exists(p$log.file)) file.remove(p$log.file)
  cat("fastsimcoal running...\n")
  err <- system2(exec, args, stdout = p$log.file)
  # err <- if(.Platform$OS.type == "unix") {
  #   system2(exec, args, stdout = p$log.file)
  # } else {
  #   shell(paste(exec, args), intern = F)
  # }
  if(err == 0) {
    p$args <- args
  } else {
    stop(
      "fastsimcoal exited with error ", err, "\n",
      "The command was:\n",
      exec, args
    )
  }
  
  p$arp.files <- NULL
  arp.files <- dir(p$label, pattern = ".arp$", full.names = TRUE)
  if(length(arp.files) > 0) {
    p$arp.files <- stats::setNames(
      arp.files, 
      paste0("rep", 1:length(arp.files))
    )
  }
  
  if(read.output) {
    if(!p$is.tpl & !is.null(p$arp.files)) {
      p <- fscRead(p)
    } else if(p$is.tpl) {
      p <- fscReadParamEst(p)
    }
  }
  if(delete.files) fscCleanup(p)
  options(opt)
  
  invisible(p)
}


# ---- Reading ----

#' @noRd
.fscMapArpLocusInfo <- function(p) {
  # expand genetic info in input parameters to matrix
  locus.info <- do.call(rbind, p$genetics)
  if(!attr(p$genetics, "chrom.diff")) {
    locus.info <- do.call(
      rbind, 
      lapply(1:attr(p$genetics, "num.chrom"), function(i) {
        cbind(chromosome = i, locus.info)
      })
    )
  }
  
  # create list to associate rows in locus info to columns read from .arp file
  mat.col <- lapply(1:nrow(locus.info), function(i) NA)
  mat.col[[1]] <- 3
  for(i in 2:nrow(locus.info)) {
    # repeat column number for consecutive DNA entries
    if(locus.info$fsc.type[i] == "DNA" & locus.info$fsc.type[i - 1] == "DNA") {
      mat.col[[i]] <- mat.col[[i - 1]]
    } else {
      next.col <- max(mat.col[[i - 1]]) + 1
      mat.col[[i]] <- if(locus.info$fsc.type[i] == "DNA") { # DNA takes 1 column
        next.col
      } else { # MICROSAT and STANDARD take num.markers columns
        next.col:(next.col + locus.info$num.markers[i] - 1)
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
        "C", swfscMisc::zero.pad(.data$chromosome),
        "_B", swfscMisc::zero.pad(.data$block), 
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
  locus.info %>% 
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
    dplyr::ungroup()
}


#' @rdname fastsimcoal
#' @export
#' 
fscReadArpFile <- function(file) {
  # read .arp file
  f <- readLines(file)
  
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


#' @rdname fastsimcoal
#' @export
#' 
fscRead <- function(p) {
  locus.info <- .fscMapArpLocusInfo(p)
  
  # loop through each .arp file
  p$genetic.data <- sapply(p$arp.files, function(file) {
    data.mat <- fscReadArpFile(file)
    
    # reform genetic data to map linkage blocks to each column
    gen.data <- lapply(1:nrow(locus.info), function(i) {
      cols <- data.mat[, locus.info$mat.col[[i]]]
      # extract DNA sequence from string in column
      if(locus.info$fsc.type[i] == "DNA") {
        cols <- substr(cols, locus.info$dna.start[i], locus.info$dna.end[i])
        if(locus.info$actual.type[i] == "SNP") {
          cols <- do.call(rbind, strsplit(cols, ""))
        }
      }
      cols <- cbind(cols)
      # give suffixes to block names with more than one locus
      colnames(cols) <- if(ncol(cols) == 1) {
        locus.info$name[i] 
      } else {
        paste0(locus.info$name[i], "_L", swfscMisc::zero.pad(1:ncol(cols)))
      }
      cols
    })
    
    # create map of column numbers in gen.data for each row of locus.info
    last.col <- 2
    locus.cols <- list()
    for(mat in gen.data) {
      col.nums <- (last.col + 1):(last.col + ncol(mat))
      locus.cols <- c(locus.cols, list(col.nums))
      last.col <- max(col.nums)
    }
    names(locus.cols) <- locus.info$name
  
    gen.data <- cbind(data.mat[, 1:2], do.call(cbind, gen.data))
    attr(gen.data, "locus.cols") <- locus.cols
    gen.data
  }, simplify = FALSE, USE.NAMES = TRUE)
  
  p$locus.info <- locus.info %>% 
    dplyr::select(.data$name, .data$chromosome, .data$block, dplyr::everything()) %>% 
    dplyr::select(-.data$mat.col, -.data$dna.col, -.data$dna.end, -.data$dna.start) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  
  invisible(p)
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
fscReadParamEst <- function(p) {
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
  
  brent.file <- readLines(brent.file) 
  brent.file <- brent.file[grep("^Param|^[[:digit:]]", brent.file)]
  brent.file <- strsplit(brent.file, "\t")
  brent <- as.data.frame(do.call(rbind, lapply(brent.file[-1], as.numeric)))
  colnames(brent) <- brent.file[[1]][1:ncol(brent)]
  
  sfs.vec <- if(length(sfs.file) == 1) {
    fscReadVector(sfs.file)
  } else {
    sapply(sfs.file, fscReadVector, simplify = FALSE)
  }
  
  p$est.params <- list(
    ml.params = fscReadVector(best.file),
    brent.lhoods = brent,
    sfs = sfs.vec
  )
  
  invisible(p)
}

#' @rdname fastsimcoal
#' @export
#' 
fscExtractLoci <- function(p, type = "all", sep.chrom = FALSE, chrom = NULL) {
  if(is.null(p$locus.info) | is.null(p$genetic.data)) {
    p <- fscRead(p)
  }
  
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
  
  .extractLocCols <- function(loc.names, gen.mat) {
    loc.cols <- unlist(attr(gen.mat, "locus.cols")[loc.names]) 
    gen.mat[, c(1:2, loc.cols), drop = FALSE]
  }
  
  # extract columns for selected chromosomes and marker types
  sapply(p$genetic.data, function(gen.mat) {
    if(sep.chrom) { # extract for each chromosome
      locus.info$chrom.label <- regmatches(
        locus.info$chromosome,
        regexpr("^C[[:alnum:]]+", locus.info$chromosome)
      )
      tapply( 
        locus.info$name, 
        locus.info$chrom.label,
        .extractLocCols,
        gen.mat = gen.mat
      )
    } else { # extract selected loci from data frame
      .extractLocCols(locus.info$name, gen.mat)
    }
  }, simplify = FALSE, USE.NAMES = TRUE)
}


#' @rdname fastsimcoal
#' @export
#' 
fscCleanup <- function(p) {
  unlink(p$label, recursive = TRUE, force = TRUE)
  if(file.exists(p$in.file)) file.remove(p$in.file)
  if(!is.null(p$est.file)) if(file.exists(p$est.file)) file.remove(p$est.file)
  if(file.exists("seed.txt")) file.remove("seed.txt")
  invisible()
}


# ---- gtypes ----

#' @rdname fastsimcoal
#' @export
#' 
fsc2gtypes <- function(p, type = c("snp", "microsat", "standard", "dna"),
                       chrom = NULL, drop.mono = TRUE) {
  # check arguments
  type <- match.arg(type)
  ploidy <- attr(p$genetics, "ploidy")
  if(ploidy == 1 & type != "dna") {
    stop("Can't create a haploid `gtypes` object for `type` other than 'dna'.")
  }
  sep.chrom <- ifelse(type == "dna" & ploidy == 1, TRUE, FALSE)
  
  # extract requested data
  gen.data <- fscExtractLoci(p, type, sep.chrom, chrom)
  if(drop.mono & ploidy > 1) {
    for(i in 1:length(gen.data)) {
      not.mono <- apply(gen.data[[i]][, -(1:2), drop = FALSE], 2, function(loc) {
        length(unique(loc)) > 1
      })
      gen.data[[i]] <- gen.data[[i]][, c(1:2, which(not.mono) + 2), drop = FALSE]
    }
  }
      
  # create list of gtypes for each replicate
  fsc.g <- purrr::imap(gen.data, function(rep.mat, rep.num) {
    if(ncol(rep.mat) <= 2) return(NULL) 
    description <- paste(p$label, rep.num, sep = ".")
    if(ploidy > 1) { # compile non-haploid genotypes
      rep.mat %>% 
        dplyr::as_tibble() %>% 
        dplyr::mutate(ind = rep(1:(dplyr::n() / ploidy), each = ploidy)) %>% 
        tidyr::gather("locus", "allele", -.data$ind, -.data$id, -.data$deme) %>% 
        dplyr::group_by(.data$deme, .data$ind, .data$locus) %>% 
        dplyr::mutate(
          allele.num = 1:dplyr::n(),
          id = paste(.data$id, collapse = "|")
        ) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(
          locus = paste(.data$locus, .data$allele.num, sep = ".")
        ) %>% 
        dplyr::select(-.data$allele.num, -.data$ind) %>% 
        tidyr::spread(.data$locus, .data$allele) %>% 
        df2gtypes(ploidy = ploidy, description = description)
    } else { # create list of sequences
      dna.list <- sapply(rep.mat, function(chrom.mat) {
        apply(
          as.matrix(chrom.mat[, -(1:2)]),
          1,
          paste,
          collapse = ""
        ) %>% 
          strsplit(split = "") %>% 
          stats::setNames(chrom.mat[, "id"]) %>% 
          ape::as.DNAbin
      }, simplify = FALSE, USE.NAMES = TRUE)
      rep.mat <- rep.mat[[1]][, c(1, 2, rep(1, length(dna.list)))]
      colnames(rep.mat)[3:ncol(rep.mat)] <- names(dna.list)
      df2gtypes(
        rep.mat, 
        ploidy = 1, 
        sequences = dna.list, 
        description = description
      ) %>% 
        labelHaplotypes
    }
  })
  
  some.null <- any(sapply(fsc.g, is.null))
  if(some.null) {
    warning("Some replicates had no specified loci - gtypes for those are NULL.")
  }
  
  fsc.g
}