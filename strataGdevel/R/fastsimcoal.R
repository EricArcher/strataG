#' @name fastsimcoal
#' @title Run fastsimcoal
#' @description Run a fastsimcoal simualtion and load results into a 
#'   \linkS4class{gtypes} object.
#'
#' @param num.pops number of populations.
#' @param Ne effective population size.
#' @param sample.size number of samples to take.
#' @param sample.time time to draw samples.
#' @param growth.rate growth rate of populations.
#' @param mig.rates list of migration matrices.
#' @param hist.ev matrix of historical events.
#' @param locus.params a list of locus parameter matrices with one element per 
#'   chromosome.
#' @param num.chrom number of chromosomes if multiple chromosomes with the same 
#'   structure are desired. Leave as \code{NULL} to take the number from the 
#'   number of elements in \code{locus.params}. If specified, the matrix in the 
#'   first element in \code{locus.params} will be used.
#' @param label character string to label files with.
#' @param num.sims number of simulations to run.
#' @param inf.site.model logical. Use infinite site model?
#' @param quiet logical. Run quietly?
#' @param delete.files logical. Delete files when done?
#' @param exec name of fastsimcoal executable.
#' 
#' @note Assumes that the program \code{fastsimcoal} is properly installed and 
#'   available on the command line. On PC's, this requires having it in a 
#'   folder in the PATH environmental variable. On Macs, the executable 
#'   should be installed in a folder like \code{/usr/local/bin}. The actual name of 
#'   the executable should be specified with the \code{exec} argument.
#' 
#' @return A list of \code{\link{gtypes}} objects for each simulated dataset.
#' 
#' @references Excoffier, L. and Foll, M (2011) fastsimcoal: a continuous-time 
#'   coalescent simulator of genomic diversity under arbitrarily complex 
#'   evolutionary scenarios Bioinformatics 27: 1332-1334.\cr
#'   \url{http://cmpg.unibe.ch/software/fastsimcoal2/}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
NULL


# check the locus.params list
.fsc.check.locus.params <- function(data.type, locus.params) {
  req.col <- switch(data.type, DNA = 4, MICROSAT = 5, SNP = 3, STANDARD = 3, FREQ = 3)
  if(!is.list(locus.params)) stop("'locus.params': not a list")
  for(x in locus.params) {
    if(!is.matrix(x)) stop("'locus.params': some elements are not matrices")
    # check data type row lengths
    for(i in nrow(x)) {
      num.col <- length(x[i, ])
      if(num.col != req.col) stop(
        paste("'locus.params': a", data.type, "row has", num.col, "parameters, but should have", req.col)
      )
    }
  }
}


# write input file
.fsc.write <- function(num.pops, Ne, sample.size, sample.time, growth.rate, 
                       mig.rates, hist.ev, data.type, locus.params, num.chrom, label) {
  
  # Check argument format
  if(!is.numeric(hist.ev) & is.matrix(hist.ev)) stop("'hist.ev' is not a numerical matrix")
  .fsc.check.locus.params(data.type, locus.params)
  
  file <- paste(label, ".par", sep = "")
  write(paste("//  <<", label, ">>  (input from 'strataG::fastsimcoal')"), file)
  write(paste(num.pops, "populations to sample"), file, append = T)
  
  write("//Population effective sizes", file, append = T)
  for(i in 1:length(Ne)) write(Ne[i], file, append = T)
  
  write("//Sample sizes", file, append = T)
  if(is.null(sample.size)) sample.size <- rep(Ne, length(num.pops))
  if(is.null(sample.time)) sample.time <- rep(0, length(num.pops))
  sample.size <- paste(sample.size, sample.time)
  for(i in 1:length(sample.size)) write(sample.size[i], file, append = T)
  
  write("//Growth rates", file, append = T)
  if(is.null(growth.rate)) growth.rate <- rep(0, num.pops)
  for(i in 1:length(growth.rate)) write(growth.rate[i], file, append = T)
  
  write("//Number of migration matrices", file, append = T)
  write(length(mig.rates), file, append = T)
  if(!is.null(mig.rates)) {
    for(i in 1:length(mig.rates)) {
      write("//migration matrix", file, append = T)
      for(r in 1:nrow(mig.rates[[i]])) write(mig.rates[[i]][r, ], file, append = T)
    }
  }
  
  write("//Historical events: time, source, sink, migrants, new size, growth rate, migr. matrix", file, append = T)
  write(ifelse(is.null(hist.ev), 0, nrow(hist.ev)), file, append = T)
  if(!is.null(hist.ev)) {
    for(i in 1:nrow(hist.ev)) write(paste(hist.ev[i, ], collapse = " "), file, append = T)
  }
  
  chrom.diff <- if(is.null(num.chrom)) {
    num.chrom <- length(locus.params)
    1
  } else {
    locus.params <- locus.params[1]
    0
  }
  write("//Number of independent loci [chromosome]", file, append = T)
  write(paste(num.chrom, chrom.diff), file, append = T)
  for(chrom in locus.params) {
    write("//Per chromosome: Number of linkage blocks", file, append = T)
    write(nrow(chrom), file, append = T)
    write("//Per block: data type, num loci, rec. rate and mut rate + optional parameters", file, append = T)
    for(i in 1:nrow(chrom)) {
      chrom.line <- paste(chrom[i, ], collapse = " ")
      chrom.line <- paste(data.type, chrom.line)
      write(chrom.line, file, append = T)
    }
  }
  
  file
}


.fsc.parse.simparam <- function(folder) {
  sp <- dir(folder, pattern = ".simparam", full.names = TRUE)
  sp <- readLines(sp)
  start <- grep("Number of linkage blocks", sp)
  end <- c(start - 1, length(sp))
  end <- end[-1]
  lapply(1:length(start), function(i) {
    chrom.sp <- sp[start[i]:end[i]]
    info <- grep("partially linked", chrom.sp, value = TRUE)
    info <- regmatches(info, gregexpr("[[:digit:]]", info))
    as.numeric(sapply(info, paste, collapse = ""))
  })
}


# parse dna sequence data and return gtypes
.fsc.parse.dna <- function(f, label, pop.data) {
  # get polymorphic positions
  poly.pos <- grep("polymorphic positions on chromosome", f) + 1
  num.poly <- f[poly.pos - 1]
  num.poly <- regmatches(num.poly, gregexpr("^#[[:space:]][[:digit:]]+", num.poly))
  num.poly <- sapply(num.poly, function(x) as.numeric(gsub("# ", "", x)))
  poly.pos <- gsub("#", "", f[poly.pos])
  poly.pos <- strsplit(poly.pos, ",")
  poly.pos[num.poly == 0] <- NA
  poly.pos <- lapply(poly.pos, as.numeric)
  
  seq.len <- .fsc.parse.simparam(label)
  seq.mat <- do.call(rbind, strsplit(pop.data[, 3], ""))
  
  dna.seqs <- lapply(1:length(seq.len), function(i) {
    seq.mat.i <- matrix("A", nrow = nrow(pop.data), ncol = sum(seq.len[[i]]))
    if(num.poly[i] > 0) {
      for(j in 1:length(poly.pos[[i]])) {
        seq.mat.i[, poly.pos[[i]][j]] <- seq.mat[, j]
      }
    }
    seqs <- lapply(1:nrow(seq.mat.i), function(j) tolower(seq.mat.i[j, ]))
    names(seqs) <- pop.data[, 2]
    as.DNAbin(seqs)
  })
  names(dna.seqs) <- paste("chrom", 1:length(dna.seqs), sep = ".")
  
  #     seq.i <- if(num.poly[i] == 0) {
  #       full.seq <- paste(rep("A", seq.len[i]), collapse = "")
  #       rep(full.seq, nrow(pop.data))
  #     } else { # otherwise add A's to pad out to full sequence length
  #       padding <- paste(rep("A", seq.len[i] - num.poly[i]), collapse = "")
  #       seq.i <- substr(pop.data[, 3], start[i], end[i])
  #       sapply(seq.i, function(x) paste(x, padding, sep = "", collapse = ""))
  #     }
  #     names(seq.i) <- pop.data[, 2]
  #     as.DNAbin(strsplit(tolower(seq.i), ""))
  #   })
  
  dna.seqs <- new("multidna", dna.seqs)
  sequence2gtypes(dna.seqs, strata = pop.data[, 1], description = label)
}


# parse microsatellite or snp data and return gtypes
.fsc.parse.msats.snps <- function(f, pop.data){
  # get diploid data
  n.loc <- ncol(pop.data) - 2
  pop.data <- do.call(
    rbind, lapply(seq(1, nrow(pop.data), 2), function(i) {
      ind <- pop.data[c(i, i + 1), ]
      locus.data <- as.vector(ind[, -(1:2)])
      c(ind[1, 1], paste(ind[, 2], collapse = "/"), locus.data)
    })
  )
  locus.data <- pop.data[, -c(1:2)]
  collapsed.loci <- do.call(
    cbind, lapply(seq(2, ncol(locus.data), by = 2), function(i) {
      a1 <- locus.data[, i - 1]
      a2 <- locus.data[, i]
      paste(a1, a2, sep = "/")
    })
  )
  colnames(collapsed.loci) <- paste("Locus", 1:ncol(collapsed.loci), sep = ".")
  rownames(collapsed.loci) <- pop.data[, 2]
  g <- df2gtypes(pop.data, ploidy = 2, id.col = 2, strata.col = 1, description = label)
  labelHaplotypes(g)$gtypes
}


# read arlequin output
.fsc.read <- function(arl.files, data.type, label) {  
  fs.gtypes <- lapply(arl.files, function(file) {
    f <- readLines(file)
    
    # get start and end points of data blocks
    start <- grep("SampleData=", f) + 1
    end <- which(f == "}") - 2
    pos <- cbind(start, end)
    
    # compile data for each population
    pop.data <- do.call(rbind, lapply(1:nrow(pos), function(i) {
      f.line <- f[pos[i, 1]:pos[i, 2]]
      f.line <- gsub("[[:space:]]+", "--", f.line)
      data.mat <- do.call(rbind, strsplit(f.line, "--"))[, -2, drop = FALSE]
      data.mat <- cbind(rep(paste("Sample", i), nrow(data.mat)), data.mat)
    }))
    
    # get data type and parse data
    data.type <- gsub("\tDataType=", "", grep("DataType=", f, value = TRUE))
    is.seq <- c(DNA = TRUE, MICROSAT = FALSE, STANDARD = FALSE)[data.type]
    if(is.seq) {
      .fsc.parse.dna(f, label, pop.data) 
    } else {
      .fsc.parse.msats.snps(f, pop.data)
    }
  })
  names(fs.gtypes) <- basename(arl.files)
  fs.gtypes
}


#' @rdname fastsimcoal
#' @export
#' 
fastsimcoal <- function(num.pops, Ne, sample.size = NULL, 
                        sample.time = NULL, growth.rate = NULL, mig.rates = NULL, hist.ev = NULL, 
                        data.type = c("DNA", "MICROSAT", "SNP", "STANDARD", "FREQ"), locus.params = NULL, 
                        num.chrom = NULL, label = "strataG_fastsimcoal", num.sims = 1, 
                        inf.site.model = TRUE, quiet = TRUE, delete.files = TRUE, exec = "fsc252") {
  
  data.type <- match.arg(data.type)
  
  # Write input file
  infile <- .fsc.write(
    num.pops = num.pops, Ne = Ne, sample.size = sample.size, sample.time = sample.time, 
    growth.rate = growth.rate, mig.rates = mig.rates, hist.ev = hist.ev, 
    data.type = data.type, locus.params = locus.params, num.chrom = num.chrom, 
    label = label
  ) 
  
  # Check/setup folder structure
  if(file.exists(label)) for(f in dir(label, full.names = T)) file.remove(f)
  
  # Run fastsimcoal
  cmd <- paste(
    exec, " -i", infile, "-n", num.sims, 
    ifelse(inf.site.model, "-I", ""), ifelse(quiet, "-q", "")
  )
  err <- system(cmd, intern = F)
  if(err != 0) stop(paste("fastsimcoal returned error code", err)) 
  arl.files <- dir(label, pattern = ".arp", full.names = T)
  if(length(arl.files) == 0) stop("fastsimcoal did not generate output")
  
  # Read and parse output
  fsc.gtypes <- .fsc.read(arl.files, data.type, label)
  
  # Cleanup
  if(delete.files) {
    unlink(label, recursive = TRUE, force = TRUE)
    file.remove(infile)
  }
  
  fsc.gtypes
}