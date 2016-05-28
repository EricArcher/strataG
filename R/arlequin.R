#' @name arlequin
#' @title Read and Write Arlequin Files
#' @description Read and write an Arlequin-formatted files.
#' 
#' @param file filename for output file.
#' @param g a \linkS4class{gtypes} object.
#' @param label label for filename(s). Default is the gtypes description if present.
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
#' @importFrom reshape2 dcast
#' @export
#' 
read.arlequin <- function(file) {
  # arp <- scan(file, what = "character", quiet = TRUE)
  # hd <- grep("HaplotypeDefinition", arp, ignore.case = TRUE)[1]
  # left.brace <- grep("[{]", arp)
  # right.brace <- grep("[}]", arp)
  # seq.start <- min(left.brace[left.brace >= hd]) + 1
  # seq.end <- min(right.brace) - 1
  # seq.lt <- strsplit(arp[seq(seq.start + 1, seq.end, 2)], "")
  # names(seq.lt) <- toupper(arp[seq(seq.start, seq.end - 1, 2)])
  # seq.lt <- lapply(seq.lt, tolower)
  # seq.lt <- lapply(seq.lt, function(x) {
  #   x[x == "?"] <- "n"
  #   x
  # })
  # sample.name <- grep("SampleName", arp, ignore.case = TRUE)
  # left.brace <- left.brace[left.brace > seq.end + 1]
  # right.brace <- right.brace[right.brace > seq.end + 1]
  # eq.symbol <- grep("=", arp)
  # freq.df <- do.call(rbind, lapply(sample.name, function(i) {
  #   strata.start <- min(eq.symbol[eq.symbol >= i])
  #   strata.end <- min(eq.symbol[eq.symbol > strata.start]) - 1
  #   strata <- paste(arp[strata.start:strata.end], collapse = " ")
  #   strata <- gsub("^[[:alnum:]=[:blank:]]*[\"]", "", strata)
  #   strata <- gsub("[\"]", "", strata)
  #   freq.start <- min(left.brace[left.brace > i]) + 1
  #   freq.end <- min(right.brace[right.brace > i]) - 1
  #   haps <- arp[seq(freq.start, freq.end - 1, by = 2)]
  #   freqs <- as.numeric(arp[seq(freq.start + 1, freq.end, by = 2)])
  #   data.frame(strata = strata, haps = haps, freqs = freqs, 
  #              stringsAsFactors = FALSE)
  # }))
  # freq.df <- dcast(freq.df, haps ~ strata, sum, value.var = "freqs")
  # freq.df$haps <- toupper(freq.df$haps)
  # list(freq.df = freq.df, dna.seq = seq.lt)
}

.writeArlequinHeader <- function(g, file) {
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
  data.type
}

.writeArlequinSequences <- function(st, file, locus) {  
  hap.freqs <- alleleFreqs(st)[[locus]][, "freq"]
  hap.freqs <- cbind(names(hap.freqs), hap.freqs)
  dna <- NULL
  if(!is.null(sequences(st))) {
    dna <- toupper(as.character(as.matrix(getSequences(sequences(st), locus))))
    dna <- apply(dna, 1, paste, collapse = "")[rownames(hap.freqs)]
  }
  hap.mat <- cbind(hap.freqs, dna)
  write(t(hap.mat), ncolumns = ncol(hap.mat), sep = "\t", file = file, append = TRUE)
}

.writeArlequinMsats <- function(st, file) {
  loc.mat <- do.call(rbind, lapply(indNames(st), function(id) {
    id.mat <- sapply(loci(st, ids = id), function(x) {
      x <- as.character(x)
      x[is.na(x)] <- "?"
      x
    })
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

.writeArlequinStructure <- function(g, file, data.type) {
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
}

#' @rdname arlequin
#' @export
#' 
write.arlequin <- function(g, label = NULL, locus = 1) {
  label <- .getFileLabel(g, label)
  file <- paste(label, ".arp", sep = "")
  
  
  data.type <- .writeArlequinHeader(g, file)
  write("[Data]", file = file, append = TRUE)
  write("[[Samples]]", file = file, append = TRUE)  
  for(st in strataSplit(g)) {
    write(paste("SampleName=\"", strata(st)[1], "\"", sep = ""), file = file, append = TRUE)
    write(paste("SampleSize=", nInd(st), sep = ""), file = file, append = TRUE)
    write("SampleData={", file = file, append = TRUE)
    if(ploidy(g) == 1) {
      .writeArlequinSequences(st, file, locus)
    } else {
      .writeArlequinMsats(st, file)
    }    
    write("}", file = file, append = TRUE)
  }
  .writeArlequinStructure(g, file, data.type)
  
  invisible(NULL)
}