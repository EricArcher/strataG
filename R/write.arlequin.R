#' @title Write Arlequin Input Files
#' @description Create an Arlequin-formatted input file from a \code{gtypes} object.
#' 
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
#' @name write.arlequin
#' @aliases arlequin
#' @importFrom reshape2 dcast
#' 
NULL

#' @rdname write.arlequin
#' @keywords internal
#' 
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

#' @rdname write.arlequin
#' @keywords internal
#' 
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

#' @rdname write.arlequin
#' @keywords internal
#' 
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

#' @rdname write.arlequin
#' @keywords internal
#' 
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

#' @rdname write.arlequin
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