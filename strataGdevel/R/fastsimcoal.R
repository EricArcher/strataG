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
#' @param hist.ev historical events.
#' @param num.chrom number of chromosomes.
#' @param locus.params locus parameters.
#' @param label character string to label files with.
#' @param num.sims number of simulations to run.
#' @param inf.site.model logical. Infinite site model?
#' @param quiet logical. Run quietly?
#' @param delete.files logical. Delete files when done?
#' @param exec fastsimcoal executable
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
#' @export
#' 
fastsimcoal <- function(num.pops, Ne, sample.size = NULL, 
  sample.time = NULL, growth.rate = NULL, mig.rates = NULL, hist.ev = NULL, 
  num.chrom = 1, locus.params = NULL, 
  label = "strataG_fastsimcoal", num.sims = 1, 
  inf.site.model = TRUE, quiet = TRUE, delete.files = TRUE,
  exec = "fsc252") {
  
  # Write input file
  hist.ev <- if(is.list(hist.ev)) do.call(rbind, hist.ev) else rbind(hist.ev)
  locus.params <- if(is.list(locus.params)) {
    do.call(rbind, locus.params) 
  } else {
    rbind(locus.params)
  }
  if(nrow(locus.params) == 1 & num.chrom > 1) {
    locus.params <- do.call(rbind, lapply(1:num.chrom, function(i) {
      locus.params[1, ]
     }))
  }
  
  file <- paste(label, ".par", sep = "")
  write(paste("//  <<", label, 
              ">>  (input from 'fastsimcoal.skeleSim.run')"), file)
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
      for(r in 1:nrow(mig.rates[[i]])) {
        write(mig.rates[[i]][r, ], file, append = T)
      }
    }
  }
  
  write("//Historical events: time, source, sink, migrants, new size, growth rate, migr. matrix", file, append = T)
  write(ifelse(is.null(hist.ev), 0, nrow(hist.ev)), file, append = T)
  if(!is.null(hist.ev)) {
    for(i in 1:nrow(hist.ev)) {
      write(paste(hist.ev[i, ], collapse = " "), file, append = T)
    }
  }
  
  write("//Number of independent loci [chromosome]", file, append = T)
  write(paste(num.chrom, "0"), file, append = T)
  for(i in 1:num.chrom) {
    write("//Per chromosome: Number of linkage blocks", file, append = T)
    write("1", file, append = T)
    write("//Per block: data type, num loci, rec. rate and mut rate + optional parameters", file, append = T)
    write(paste(locus.params[i, ], collapse = " "), file, append = T)
  }
  
  # Check/setup folder structure
  if(file.exists(label)) for(f in dir(label, full.names = T)) file.remove(f)
  
  # Run fastsimcoal
  cmd <- paste(exec, " -i", file, "-n 1",
               ifelse(inf.site.model, "-I", ""), ifelse(quiet, "-q", "")
  )
  err <- system(cmd, intern = F)
  if(err != 0) stop(paste("fastsimcoal returned error code", err)) 
  arl.files <- dir(label, pattern = ".arp", full.names = T)
  if(length(arl.files) == 0) stop("fastsimcoal did not generate output")
  
  # Read arlequin output
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
      data.mat <- do.call(rbind, strsplit(f.line, "--"))[, -2]
      data.mat <- cbind(rep(paste("Sample", i), nrow(data.mat)), data.mat)
    }))
    
    # get data type
    data.type <- f[grep("DataType=", f)]
    data.type <- gsub("\tDataType=", "", data.type)
    is.seq <- switch(data.type, DNA = T, MICROSAT = F, STANDARD = F)
    if(is.seq) {    
      # get sequence length markers
      poly.pos.lines <- grep("polymorphic positions on chromosome", f, value = T)
      num.poly <- as.numeric(sapply(strsplit(poly.pos.lines, " "), function(x) x[2]))
      end <- cumsum(num.poly)
      start <- c(1, end[1:(length(end) - 1)] + 1)
      
      seq.len <- locus.params[, 1]
      dna.seqs <- lapply(1:length(seq.len), function(i) {
        seq.i <- if(num.poly[i] == 0) {
          full.seq <- paste(rep("A", seq.len[i]), collapse = "")
          rep(full.seq, nrow(pop.data))
        } else { # otherwise add A's to pad out to full sequence length
          padding <- paste(rep("A", seq.len[i] - num.poly[i]), collapse = "")
          seq.i <- substr(pop.data[, 3], start[i], end[i])
          sapply(seq.i, function(x) paste(x, padding, sep = "", collapse = ""))
        }
        names(seq.i) <- pop.data[, 2]
        as.DNAbin(strsplit(tolower(seq.i), ""))
      })
      dna.seqs <- new("multidna", dna.seqs)
      sequence2gtypes(dna.seqs, strata = pop.data[, 1], description = label)
      
#       # replace sequence with all A's if there are no variable sites
#       n.loc <- locus.params[1, 1]
#       if(pop.data[1, 3] == "?") {
#         full.seq <- paste(rep("A", n.loc), collapse = "")
#         pop.data[, 3] <- rep(full.seq, nrow(pop.data))
#       } else { # otherwise add A's to pad out to full sequence length
#         partial.seq <- paste(rep("A", n.loc - nchar(pop.data[1, 3])), 
#                              collapse = "")
#         pop.data[, 3] <- sapply(pop.data[, 3], function(x) {
#           paste(x, partial.seq, sep = "", collapse = "")
#         })
#       }
#       dna.seq <- strsplit(pop.data[, 3], "")
#       names(dna.seq) <- pop.data[, 2]
#       sequence2gtypes(dna.seq, strata = pop.data[, 1], description = label)
    } else {
      # compile diploid data
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
      df2gtypes(pop.data, ploidy = 2, id.col = 2, strata.col = 1, description = label)
    }
  })
  names(fs.gtypes) <- basename(arl.files)
  
  if(delete.files) file.remove(c(dir(label, full.names = T), label, file))
  
  fs.gtypes
}
