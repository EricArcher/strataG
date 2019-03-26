#' @title Run fastsimcoal
#' @description Run a fastsimcoal simulation and load results into a 
#'   \linkS4class{gtypes} object.
#'
#' @param label character string to label files with.
#' @param seed random number seed for simulation.
#' @param exec name of fastsimcoal executable.
#' @param num.cores number of cores to use.
#' @param p list of fastsimcoal input parameters and output.
#' @param num.sims number of simulation replicates to run.
#' @param dna.to.snp convert DNA sequences to numerical SNPs?
#' @param max.snps maximum number of SNPs to retain.
#' @param all.sites retain all sites? If \code{FALSE}, only polymorphic DNA 
#'   sites will be returned. This includes SNP blocks as they are simulated as 
#'   DNA sequences.
#' @param inf.sites use infinite sites model? If \code{TRUE}, all mutations are 
#'   retained in the output, thus the number of sites for SNPs or DNA sequences 
#'   will potentially be greater than what was requested.
#' @param no.arl.output do not output arlequin files.
#' @param trees output NEXUS formatted coalescent trees for all replicates?
#' @param sfs.type type of site frequency spectrum to compute for each 
#'   population sample: `daf` = derived allele frequency (unfolded), 
#'   `maf` = minor allele frequency (folded).
#' @param jobs output one simulated or bootstrapped SFS per file in a 
#'   separate directory?
#' @param nonpar.boot number of bootstraps to perform on polymorphic sites to
#'   extract SFS.
#' @param maxlhood Perform parameter estimation by maximum composite 
#'   likelihood from the SFS?
#' @param num.loops number of loops (ECM cycles) to be performed when 
#'   estimating parameters from SFS. Default is 20.
#' @param min.num.loops number of loops (ECM cycles) for which the 
#'   likelihood is computed on both monomorphic and polymorphic sites. Default 
#'   is 20.
#' @param brent Tolerance level for Brent optimization.   
#'   Smaller value imply more precise estimations, but require more 
#'   computation time. Default = 0.01. Value is restricted between 
#'   1e-5 and 1e-1.
#' @param save.est do not delete .est parameter estimation files during cleanup?
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
#' @name fscRun
#' @export
#' 
fscRun <- function(p, num.sims = 1, dna.to.snp = FALSE, max.snps = NULL, 
                   all.sites = TRUE, inf.sites = FALSE, 
                   sfs.type = c("maf", "daf"), jobs = FALSE, nonpar.boot = NULL, 
                   no.arl.output = FALSE, maxlhood = FALSE, num.loops = 20, 
                   min.num.loops = 20, brent = 0.01, trees = FALSE,
                   num.cores = NULL, seed = NULL, exec = "fsc26") {
  
  run.params <- as.list(environment())
  run.params$p <- NULL
  
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
  
  args <- c(
    ifelse(p$is.tpl, "--tplfile", "--ifile"), p$files$in.file,
    "--numsims", num.sims
  )
  if(p$is.tpl) args <- c(args, "-e", p$files$est.file)
  if(!is.null(seed)) args <- c(args, "--seed ", seed)
  if(trees) args <- c(args, "--tree")
  if(dna.to.snp) {
    args <- c(args, "--dnatosnp", ifelse(is.null(max.snps), 0, max.snps))
    args <- c(args, switch(match.arg(sfs.type), maf = "--msfs", daf = "--dsfs"))
  }
  if(jobs) args <- c(args, "--jobs")
  if(!is.null(nonpar.boot)) c(args, "--nonparboot", nonpar.boot)
  if(all.sites) args <- c(args, "--allsites")
  if(inf.sites) args <- c(args, "--inf")
  if(no.arl.output) args <- c(args, "--noarloutput")
  if(p$is.tpl) {
    args <- c(args, "--maxlhood")
    args <- c(args, "--numloops", num.loops)
    args <- c(args, "--minnumloops", min.num.loops)
    args <- c(
      args, 
      "--brentol", 
      ifelse(brent > 1e-1, 1e-1, ifelse(brent < 1e-5, 1e-5, brent))
    )
  }
  args <- c(args, cores.spec)
  p$run.params <- run.params
  p$run.params$args <- paste(args, collapse = " ")
  
  p$files$log.file <- file.path(p$label, paste0(p$label, ".log"))
  
  if(!dir.exists(p$label)) dir.create(p$label)
  if(file.exists(p$files$log.file)) file.remove(p$files$log.file)
  cat(format(Sys.time()), "running fastsimcoal...\n")
  err <- system2(exec, p$run.params$args, stdout = p$files$log.file)
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
  p$files$arp.files <- NULL
  arp.files <- dir(p$label, pattern = ".arp$", full.names = TRUE)
  if(length(arp.files) > 0) {
    arp.files <- arp.files[order(nchar(arp.files), arp.files)]
    p$files$arp.files <- stats::setNames(
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


#' @noRd
#' 
.fscMapArpLocusInfo <- function(p) {
  # expand genetic info in input parameters to matrix
  locus.info <- do.call(rbind, p$settings$genetics)
  if(!attr(p$settings$genetics, "chrom.diff")) {
    num.blocks <- nrow(locus.info)
    num.chrom <- attr(p$settings$genetics, "num.chrom")
    locus.info <- replicate(num.chrom, locus.info, simplify = FALSE)
    locus.info <- do.call(rbind, locus.info)
    locus.info <- cbind(
      chromosome = rep(1:num.chrom, each = num.blocks),
      locus.info
    )
  }

  # associate rows in locus info to columns in .arp file (if all.sites = TRUE)
  prev.type <- dplyr::lag(locus.info$fsc.type)
  new.col <- as.numeric(!(locus.info$fsc.type == "DNA" & prev.type == "DNA"))
  new.col <- new.col * ifelse(
    locus.info$fsc.type == "DNA", 1, locus.info$num.markers
  )
  if(locus.info$fsc.type[1] == "DNA") new.col[1] <- 1
  mat.col.end <- cumsum(new.col) + 2
  mat.col.start <- mat.col.end - ifelse(
    locus.info$fsc.type == "DNA", 0, locus.info$num.markers - 1
  )
  
  p$locus.info <- locus.info %>% 
    dplyr::group_by(.data$chromosome) %>% 
    dplyr::mutate( # by chromosome information
      block = 1:dplyr::n(),
      chrom.pos.end = cumsum(.data$num.markers),
      chrom.pos.start = .data$chrom.pos.end - .data$num.markers + 1
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      name = paste0( # create block names
        "C", .zeroPad(.data$chromosome),
        "B", .zeroPad(.data$block), 
        "_", .data$actual.type
      ),
      mat.col.start = mat.col.start,
      mat.col.end = mat.col.end
    ) %>% 
    dplyr::group_by(.data$mat.col.start) %>% 
    dplyr::mutate( # character positions of DNA loci
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
    dplyr::select(
      .data$name, .data$chromosome, .data$block, 
      .data$actual.type:.data$param.6,
      .data$chrom.pos.start, .data$chrom.pos.end,
      .data$mat.col.start, .data$mat.col.end, 
      .data$dna.start, .data$dna.end
    ) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  
  if(
    p$run.params$inf.sites & 
    p$run.params$all.sites & 
    any(p$locus.info$fsc.type == "DNA")
  ) {
    max.chrom <- max(p$locus.info$chromosome)
    max.block <- max(p$locus.info$block)
    p$locus.info <- do.call(
      rbind, 
      lapply(split(p$locus.info, p$locus.info$mat.col.start), function(data.col) {
        if(unique(data.col$fsc.type) != "DNA") return(data.col)
        data.col %>% 
          dplyr::slice(1) %>% 
          dplyr::mutate(
            name = paste0(
              "C", .zeroPad(min(data.col$chromosome), max.chrom), 
              "B", .zeroPad(min(data.col$block), max.block),
              "C", .zeroPad(max(data.col$chromosome), max.chrom),
              "B", .zeroPad(max(data.col$block), max.block), "_",
              .data$actual.type, collapse = ""
            ),
            chromosome = NA,
            block = NA,
            num.markers = sum(data.col$num.markers),
            recomb.rate = mean(data.col$recomb.rate),
            mut.rate = mean(data.col$mut.rate),
            param.5 = mean(data.col$param.5),
            chrom.pos.start = NA,
            chrom.pos.end = NA,
            dna.start = NA,
            dna.end = NA
          )
      })
    )
  }
  
  invisible(p)
}


#' @rdname fscRun
#' @export
#' 
fscCleanup <- function(label, save.est = TRUE) {
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