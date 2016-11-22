#' @title Read and Write Arlequin Files
#' @description Read an Arlequin-formatted project input file into a 
#'   \code{gtypes} object, or write an input file from a \code{gtypes} object.
#' 
#' @param file filename of an input arlequin project (.arp) file.
#'   \code{read.arlequin} can only read files with DataType of FREQUENCY, DNA, 
#'   or MICROSAT.
#' @param g a \linkS4class{gtypes} object.
#' @param label label for Arlequin project filename (.arp). If \code{NULL}, the 
#'   gtypes description is used if present.
#' @param locus numeric or character designation of which locus to write for 
#'   haploid data.
#' 
#' @references Excoffier, L.G. Laval, and S. Schneider (2005) 
#'   Arlequin ver. 3.0: An integrated software package for population genetics 
#'   data analysis. Evolutionary Bioinformatics Online 1:47-50.\cr
#'   Available at \url{http://cmpg.unibe.ch/software/arlequin3/}
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @name arlequin
#' @aliases arlequin
#' 
NULL

#' @rdname arlequin
#' @export
#' 
read.arlequin <- function(file) {
  arp <- scan(file, what = "character", sep = "\n", quiet = TRUE)
  
  # function to get value folowing "=" for a header title
  getValues <- function(title, arp) {
    x <- grep(title, arp, ignore.case = TRUE, value = TRUE)
    gsub("^[[:alnum:]]+=", "", x)
  }
  
  # get header info
  title <- gsub("[[:punct:]]+", "", getValues("Title", arp))
  data.type <- gsub("[[:punct:]]+", "", getValues("DataType", arp))
  
  if(!toupper(data.type) %in% c("FREQUNENCY", "DNA", "MICROSAT")) {
    stop("DataType must be of FREQUENCY, DNA, or MICROSAT")
  }
  
  missing.data <- gsub("[']|\"", "", getValues("MissingData", arp))
  locus.separator <- gsub("[']|\"", "", getValues("LocusSeparator", arp))
  
  # find start and stop locations for sample data
  sample.name <- gsub("[']|\"", "", getValues("SampleName", arp))
  data.start <- grep("SampleData", arp, ignore.case = TRUE) + 1
  data.end <- grep("[}]", arp)[1:length(data.start)] - 1
  
  # read genetic data into a single matrix
  split.char <- switch(
    locus.separator,
    WHITESPACE = "[[:space:]]",
    TAB = "[[:space:]]",
    NONE = "",
    locus.separator
  )
  gen.data <- do.call(rbind, lapply(1:length(data.start), function(i) {
    lines <- strsplit(arp[data.start[i]:data.end[i]], split = split.char)
    lines <- lapply(lines, function(x) {
      x <- x[x != ""]
      x[x == missing.data] <- NA
      x
    })
    max.len <- max(sapply(lines, length))
    mat <- do.call(rbind, lapply(lines, function(x) {
      if(length(x) < max.len) c(rep(NA, max.len - length(x)), x) else x
    }))
    for(r in 2:nrow(mat)) {
      if(is.na(mat[r, 1])) mat[r, 1] <- mat[r - 1, 1]
    }
    cbind(rep(sample.name[i], nrow(mat)), mat)
  }))
  
  # function to organize haplotype matrix
  haploid.mat <- function(x) {      
    new.mat <- do.call(rbind, lapply(1:nrow(x), function(i) {
      mat <- x[rep(i, as.numeric(x[i, 3])), , drop = FALSE]
      if(nrow(mat) > 1) mat[, 2] <- paste(mat[, 2], 1:nrow(mat), sep = "_")
      mat[, 2] <- paste(mat[, 2], x[i, 1], sep = "_")
      cbind(mat[, c(2, 1), drop = FALSE], rep(x[i, 2], nrow(mat)))
    }))
    colnames(new.mat) <- c("id", "strata", "gene")
    new.mat
  }
  
  # create gtypes 
  switch(
    toupper(data.type),
    "DNA" = {  
      new.mat <- haploid.mat(gen.data)
      seq.mat <- gen.data[, c(2, 4), drop = FALSE]
      seq.mat <- seq.mat[!duplicated(seq.mat[, 1]), , drop = FALSE]
      seq.mat <- seq.mat[order(seq.mat[, 1]), ]
      dna <- lapply(strsplit(seq.mat[, 2], ""), function(x) {
        x[x == missing.data] <- "n"
        tolower(x)
      })
      names(dna) <- seq.mat[, 1]
      dna <- as.DNAbin(dna)
      df2gtypes(new.mat, ploidy = 1, sequences = dna, description = title)
    },
    "FREQUENCY" = {      
      new.mat <- haploid.mat(gen.data)
      df2gtypes(new.mat, ploidy = 1, description = title)
    },
    "MICROSAT" = {
      ploidy <- unname(table(gen.data[, 2])[1])
      new.mat <- do.call(rbind, split(gen.data[, -(1:3)], gen.data[, 2]))
      freqs <- na.omit(gen.data[, 1:3])
      new.mat <- cbind(freqs, new.mat[freqs[, 2], ])
      new.mat <- do.call(rbind, lapply(1:nrow(new.mat), function(i) {
        x <- new.mat[rep(i, as.numeric(new.mat[i, 3])), , drop = FALSE]
        if(nrow(x) > 1) x[, 2] <- paste(x[, 2], 1:nrow(x), sep = "_")
        rownames(x) <- x[, 2]
        x
      }))
      new.mat <- new.mat[, c(2, 1, 4:ncol(new.mat))]
      loc.names <- paste("Locus", 1:(ncol(gen.data) - 3), sep = "")
      loc.names <- paste(rep(loc.names, each = ploidy), 1:ploidy, sep = ".")
      colnames(new.mat) <- c("id", "strata", loc.names)
      df2gtypes(new.mat, ploidy = ploidy, description = title)
    }
  )
}


#' @rdname arlequin
#' @export
#' 
write.arlequin <- function(g, label = NULL, locus = 1) {
  label <- .getFileLabel(g, label)
  file <- paste(label, ".arp", sep = "")
  
  # Header
  write("[Profile]", file = file)
  write(paste("Title=\"", description(g), "\"", sep = ""), file = file, append = TRUE)
  write(paste("NbSamples=", nStrata(g), sep = ""), file = file, append = TRUE)
  data.type <- if(ploidy(g) > 1) {
    "MICROSAT"
  } else if(is.null(sequences(g))) {
    "FREQUENCY"
  } else {
    "DNA"
  }
  write(paste("DataType=", data.type, sep = ""), file = file, append = TRUE)
  g.data <- paste("GenotypicData=", ifelse(ploidy(g) == 1, 0, 1), sep = "")
  write(g.data, file = file, append = TRUE)
  write("GameticPhase=0", file = file, append = TRUE)
  write("MissingData='?'", file = file, append = TRUE)
  write("LocusSeparator=WHITESPACE", file = file, append = TRUE)
  
  # Data
  write("[Data]", file = file, append = TRUE)
  write("[[Samples]]", file = file, append = TRUE)  
  for(st in strataSplit(g)) {
    write(paste("SampleName=\"", strata(st)[1], "\"", sep = ""), file = file, append = TRUE)
    write(paste("SampleSize=", nInd(st), sep = ""), file = file, append = TRUE)
    write("SampleData={", file = file, append = TRUE)
    if(ploidy(g) == 1) { # Sequences
      hap.freqs <- alleleFreqs(st)[[locus]][, "freq"]
      hap.freqs <- cbind(names(hap.freqs), hap.freqs)
      dna <- NULL
      if(!is.null(sequences(st))) {
        dna <- toupper(as.character(as.matrix(getSequences(sequences(st), locus))))
        dna <- apply(dna, 1, paste, collapse = "")[rownames(hap.freqs)]
      }
      hap.mat <- cbind(hap.freqs, dna)
      write(t(hap.mat), ncolumns = ncol(hap.mat), sep = "\t", file = file, append = TRUE)
    } else { # Microsats
      loc.mat <- do.call(rbind, lapply(indNames(st), function(id) {
        id.mat <- t(as.array(st, id))
        id.mat[is.na(id.mat)] <- "?"
        do.call(rbind, lapply(1:nrow(id.mat), function(i) {
          line.i <- if(i == 1) c(id, "1") else {
            id.pad <- paste(rep(" ", nchar(id)), collapse = "")
            c(id.pad, " ")
          }
          c(line.i, id.mat[i, ])
        }))
      }))
      write(t(loc.mat), ncolumns = ncol(loc.mat), sep = "\t", file = file, append = TRUE)
    }    
    write("}", file = file, append = TRUE)
  }
  
  # Structure
  write("[[Structure]]", file = file, append = TRUE)
  st.name <- paste("A group of", nStrata(g), "populations analyzed for", data.type)
  st.name <- paste("StructureName=\"", st.name, "\"", sep = "")
  write(st.name, file = file, append = TRUE)
  write("NbGroups=1", file = file, append = TRUE)
  write("Group= {", file = file, append = TRUE)
  for(st in strataNames(g)) {
    write(paste("\"", st, "\"", sep = ""), file = file, append = TRUE)
  }
  write("}", file = file, append = TRUE)
  
  invisible(NULL)
}