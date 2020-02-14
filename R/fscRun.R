#' @title Run fastsimcoal
#' @description Run a fastsimcoal simulation. 
#'
#' @param p list of fastsimcoal input parameters and output produced by 
#'   \link{fscWrite}.
#' @param num.sims number of simulation replicates to run.
#' @param dna.to.snp convert DNA sequences to numerical SNPs?
#' @param max.snps maximum number of SNPs to retain.
#' @param all.sites retain all sites? If \code{FALSE}, only polymorphic DNA 
#'   sites will be returned. This includes SNP blocks as they are simulated as 
#'   DNA sequences.
#' @param inf.sites use infinite sites model? If \code{TRUE}, all mutations are 
#'   retained in the output, thus the number of sites for SNPs or DNA sequences 
#'   will potentially be greater than what was requested.
#' @param sfs.type type of site frequency spectrum to compute for each 
#'   population sample: `daf` = derived allele frequency (unfolded), 
#'   `maf` = minor allele frequency (folded).
#' @param nonpar.boot number of bootstraps to perform on polymorphic sites to
#'   extract SFS.
#' @param no.arl.output do not output arlequin files.
#' @param num.loops number of loops (ECM cycles) to be performed when 
#'   estimating parameters from SFS. Default is 20.
#' @param min.num.loops number of loops (ECM cycles) for which the 
#'   likelihood is computed on both monomorphic and polymorphic sites. Default 
#'   is 20.
#' @param brentol Tolerance level for Brent optimization.   
#'   Smaller value imply more precise estimations, but require more 
#'   computation time. Default = 0.01. Value is restricted between 
#'   1e-5 and 1e-1.
#' @param trees output NEXUS formatted coalescent trees for all replicates?
#' @param num.cores number of cores to use. If set to \code{NULL}, the value
#'   will be what is reported by \code{\link[parallel]{detectCores} - 1}.
#' @param seed random number seed for simulation.
#' @param quiet logical indicating if fastsimcoal2 should be run in quiet mode.
#' @param exec name of fastsimcoal executable.
#' @param label character string of file run labels prefixes.
#' 
#' @return 
#' \describe{
#'  \item{fscRun}{Runs the \code{fastsimcoal2} simulation and returns a
#'    list containing run parameters and a data frame used by 
#'    \code{\link{fscRead}} to parse the genotypes generated (if 
#'    Arlequin-formatted output was requested).}
#'  \item{fscCleanup}{Deletes all files associated with the simulation 
#'    identified by \code{label}.}
#'  }
#' 
#' @note \code{fastsimcoal2} is not included with `strataG` and must be
#'   downloaded separately. Additionally, it must be installed such that it can
#'   be run from the command line in the current working directory. 
#'   The function \code{fscTutorial()} will open a detailed tutorial on the 
#'   interface in your web browser.
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
#' @seealso \code{\link{fsc.input}}, \code{\link{fscWrite}}, 
#'  \code{\link{fscRead}}
#'  
#' @examples \dontrun{
#' #' # three demes with optional names
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
#' params <- fscWrite(demes = demes, events = events, genetics = genetics)
#' 
#' # runs 100 replicates, converting all DNA sequences to 0/1 SNPs
#' # will also output the MAF site frequency spectra (SFS) for all SNP loci.
#' params <- fscRun(params, num.sim = 100, dna.to.snp = TRUE, num.cores = 3)
#' }
#' 
#' @name fscRun
#' @export
#' 
fscRun <- function(p, num.sims = 1, dna.to.snp = FALSE, max.snps = 0, 
                   sfs.type = c("maf", "daf"), nonpar.boot = NULL, 
                   all.sites = TRUE, inf.sites = FALSE, 
                   no.arl.output = FALSE, num.loops = 20, 
                   min.num.loops = 20, brentol = 0.01, trees = FALSE,
                   num.cores = 1, seed = NULL, quiet = TRUE, 
                   exec = "fsc26") {
  
  run.params <- as.list(environment())
  run.params$p <- NULL
  run.params$sfs.type <- match.arg(sfs.type)
  
  is.tpl <- !is.null(p$sim.params)
  is.est <- is.tpl & !is.null(p$files$est)
  is.def <- is.tpl & !is.null(p$files$def)
  
  if(!any(p$settings$genetics$fsc.type == "DNA") & dna.to.snp) {
    warning(
      "'dna.to.snp' set to 'FALSE' because 'fscBlock_dna()' or ",
      "'fscBlock_snp()' not used.",
      call. = FALSE
    )
    dna.to.snp <- FALSE
  }
  
  opt <- options(scipen = 999)
  
  cores.spec <- if(!is.null(num.cores)) {
    num.cores <- max(1, num.cores)
    num.cores <- min(num.cores, min(parallel::detectCores(), 12))
    if(num.cores == 1) "" else {
      paste(c("--cores", "--numBatches"), num.cores, collapse = " ")
    }
  } else ""

  args <- c(ifelse(is.tpl, "--tplfile", "--ifile"), p$files$input)
  args <- c(args, "--numsims", num.sims)
  if(is.def) {
    args <- c(args, "-f", p$files$def)
  } else if(is.est) {
    args <- c(args, "-e", p$files$est)
  }
  if(!is.null(seed)) args <- c(args, "--seed", seed)
  if(trees) args <- c(args, "--tree")
  if(dna.to.snp | is.est) {
    args <- c(args, "--dnatosnp", max.snps)
    if(!is.def) {
      args <- c(
        args, switch(match.arg(sfs.type), maf = "--msfs", daf = "--dsfs")
      )
      args <- c(args, "--jobs")
      if(!is.null(nonpar.boot)) c(args, "--nonparboot", nonpar.boot)
    }
  }
  if(all.sites) args <- c(args, "--allsites")
  if(inf.sites) args <- c(args, "--inf")
  if(no.arl.output) args <- c(args, "--noarloutput")
  if(is.est) { 
    if(brentol < 1e-5) brentol <- 1e-5
    if(brentol > 1e-1) brentol <- 1e-1
    args <- c(args, "--maxlhood")
    args <- c(args, "--numloops", num.loops)
    args <- c(args, "--minnumloops", min.num.loops)
    args <- c(args, "--brentol", brentol)
  }
  args <- c(args, cores.spec)
  if(quiet) args <- c(args, "--quiet")
  
  p$run.params <- run.params
  p$run.params$args <- paste(args, collapse = " ")
  p$files$log <- paste0(p$label, ".log")
  
  unlink(p$label, recursive = TRUE, force = TRUE)
  cat(format(Sys.time()), "running fastsimcoal2...\n")
  wd <- getwd()
  err <- tryCatch({
    setwd(p$folder)
    system2(exec, p$run.params$args, stdout = p$files$log)
  }, finally = setwd(wd))
  if(err != 0) {
    stop(
      format(Sys.time()), " fastsimcoal exited with error ", err, ". ",
      "The command was:\n", exec, " ", p$run.params$args
    )
  }
  
  # get .arp filenames and map loci
  if(!no.arl.output & is.null(p$files$est)) {
    arp.files <- NULL
    while(is.null(arp.files)) {  
      Sys.sleep(1)
      folder <- file.path(p$folder, p$label)
      arp.files <- dir(folder, pattern = ".arp$", full.names = TRUE)
      if(length(arp.files) == 0) {
        cat(format(Sys.time()), "waiting for .arp files to finish writing...\n")
        arp.files <- NULL
      }
    }
    p <- .fscMapArpLocusInfo(p)
  }

  cat(format(Sys.time()), "run complete\n")
  options(opt)
  invisible(p)
}


#' @noRd
#' 
.fscMapArpLocusInfo <- function(p) {
  # expand genetic info in input parameters to matrix
  locus.info <- p$settings$genetics
  if(!attr(p$settings$genetics, "chrom.diff")) {
    num.chrom <- attr(p$settings$genetics, "num.chrom")
    num.blocks <- nrow(locus.info)
    locus.info <- locus.info[rep(1:nrow(locus.info), num.chrom), ]
    locus.info$chromosome <- rep(1:num.chrom, each = num.blocks)
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
      lapply(
        split(p$locus.info, p$locus.info$mat.col.start), 
        function(data.col) {
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
        }
      )
    )
  }
  
  invisible(p)
}


#' @rdname fscRun
#' @param folder character string giving the root working folder where
#'   input files and output resides
#' @export
#' 
fscCleanup <- function(label, folder = ".") {
  # remove label folder
  wd <- getwd()
  tryCatch({
    setwd(folder)
    unlink(label, recursive = TRUE, force = TRUE)
    if(file.exists("seed.txt")) file.remove("seed.txt")
    files <- dir(pattern = paste0("^", label))
    if(length(files) != 0) {
      # don't remove R script files
      r.files <- grep("\\.[rR]$", files, value = TRUE)
      files <- setdiff(files, r.files)
      # remove only files, not other directories that start with label
      files <- files[utils::file_test("-f", files)]
      file.remove(files)
    }
  }, finally = setwd(wd))
  invisible(NULL)
}


#' @rdname fscRun
#' @export
#' 
fscTutorial <- function() {
  utils::browseURL(system.file("fastsimcoal2.html", package = "strataG"))
}