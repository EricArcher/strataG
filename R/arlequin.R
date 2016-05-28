#' @name arlequin
#' @title Read and Write Arlequin Files
#' @description Read and write an Arlequin-formatted files.
#' 
#' @param file filename for output file.
#' @param g a \linkS4class{gtypes} object.
#' @param label label for filename(s). Default is the gtypes description 
#'   if present.
#' @param data.type type of data. Can be "DNA", "RFLP", or "MICROSAT".
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
  arp <- scan(file, what = "character", quiet = TRUE)
  hd <- grep("HaplotypeDefinition", arp, ignore.case = TRUE)[1]
  left.brace <- grep("[{]", arp)
  right.brace <- grep("[}]", arp)
  seq.start <- min(left.brace[left.brace >= hd]) + 1
  seq.end <- min(right.brace) - 1
  seq.lt <- strsplit(arp[seq(seq.start + 1, seq.end, 2)], "")
  names(seq.lt) <- toupper(arp[seq(seq.start, seq.end - 1, 2)])
  seq.lt <- lapply(seq.lt, tolower)
  seq.lt <- lapply(seq.lt, function(x) {
    x[x == "?"] <- "n"
    x
  })
  sample.name <- grep("SampleName", arp, ignore.case = TRUE)
  left.brace <- left.brace[left.brace > seq.end + 1]
  right.brace <- right.brace[right.brace > seq.end + 1]
  eq.symbol <- grep("=", arp)
  freq.df <- do.call(rbind, lapply(sample.name, function(i) {
    strata.start <- min(eq.symbol[eq.symbol >= i])
    strata.end <- min(eq.symbol[eq.symbol > strata.start]) - 1
    strata <- paste(arp[strata.start:strata.end], collapse = " ")
    strata <- gsub("^[[:alnum:]=[:blank:]]*[\"]", "", strata)
    strata <- gsub("[\"]", "", strata)
    freq.start <- min(left.brace[left.brace > i]) + 1
    freq.end <- min(right.brace[right.brace > i]) - 1
    haps <- arp[seq(freq.start, freq.end - 1, by = 2)]
    freqs <- as.numeric(arp[seq(freq.start + 1, freq.end, by = 2)])
    data.frame(strata = strata, haps = haps, freqs = freqs, 
               stringsAsFactors = FALSE)
  }))
  freq.df <- dcast(freq.df, haps ~ strata, sum, value.var = "freqs")
  freq.df$haps <- toupper(freq.df$haps)
  list(freq.df = freq.df, dna.seq = seq.lt)
}

.writeArlequinHeader <- function(g, file, data.type) {
  write("[Profile]", file = file)
  write(paste("Title=\"", description(g), "\"", sep = ""), file = file, append = TRUE)
  write(paste("NbSamples=", nInd(g), sep = ""), file = file, append = TRUE)
  write(paste("DataType=", data.type, sep = ""), file = file, append = TRUE)
  g.data <- paste("GenotypicData=", ifelse(ploidy(g) == 1, 0, 1))
  write(g.data, file = file, append = TRUE)
  write("MissingData='?'", file = file, append = TRUE)
  write("LocusSeparator=WHITESPACE", file = file, append = TRUE)
}

.writeArlequinSequences <- function(g, file, locus) {      
  write("[[HaplotypeDefinition]]", file = file, append = TRUE)
  write("HaplListName=\"Haplotypes\"", file = file, append = TRUE)
  write("HaplList={", file = file, append = TRUE)
  dna <- as.character(as.matrix(sequences(g, locus)))
  for(x in rownames(dna)) {
    x.seq <- paste(dna[x, ], collapse = "")
    write(paste(x, x.seq), file = file, append = TRUE)
  }
  write("}", file = file, append = TRUE)
}

.writeArlequinMsats <- function(g, file) {
  write("[[Samples]]", file = file, append = TRUE)
  for(st in strataSplit(g)) {
    write(paste("SampleName=\"", strata(st)[1], "\"", sep = ""), file = file, append = TRUE)
    write(paste("SampleSize=\"", nInd(st), "\"", sep = ""), file = file, append = TRUE)
    write("SampleData={", file = file, append = TRUE)
    for(id in indNames(st)) {
      id.mat <- sapply(loci(st, ids = id), as.character)
      id.mat <- id.mat[is.na(id.mat)] <- "?"
      for(i in 1:nrow(id.mat)) {
        hdr <- if(i == 1) paste(id, "1") else ""
        als <- paste(c(hdr, id.mat[i, ]), collapse = " ")
        write(als, file = file, append = TRUE)
      }
      write("}", file = file, append = TRUE)
    }
  }
}

.writeArlequinStructure <- function(g, file, data.type) {
  write("[[Structure]]", file = file, append = TRUE)
  st.name <- paste(
    "StructureName=\"A group of", nStrata(g), "populations analyzed for", data.type
  )
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
write.arlequin <- function(g, label = NULL, data.type = c("DNA", "MICROSAT", "RFLP"),
                           locus = 1) {
  label <- .getFileLabel(g, label)
  file <- paste(label, ".arp", sep = "")
  
  data.type <- match.arg(data.type)
  if(ploidy(g) > 1 & data.type == "DNA") {
    stop("'data.type' cannot be DNA if 'g' is not haploid.")
  }
  if(ploidy(g) == 1 & data.type %in% c("RFLP", "MICROSAT")) {
    stop("'data.type' cannot be RFLP or MICROSAT if 'g' is haploid.")
  }
  
  .writeArlequinHeader(g, file, data.type)
  
  write("[Data]", file = file, append = TRUE)
  if(ploidy(g) == 1) {
    if(!is.null(sequences(g))) .writeArlequinSequences(g, file, locus)
  } else {
    .writeArlequinMsats(g, file)
  }
  .writeArlequinStructure(g, file, data.type)
  
  invisible(NULL)
}