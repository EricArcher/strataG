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
      dna <- ape::as.DNAbin(dna)
      df2gtypes(new.mat, ploidy = 1, sequences = dna, description = title)
    },
    "FREQUENCY" = {      
      new.mat <- haploid.mat(gen.data)
      df2gtypes(new.mat, ploidy = 1, description = title)
    },
    "MICROSAT" = {
      ploidy <- unname(table(gen.data[, 2])[1])
      new.mat <- do.call(rbind, split(gen.data[, -(1:3)], gen.data[, 2]))
      freqs <- stats::na.omit(gen.data[, 1:3])
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
  if(is.numeric(locus)) locus <- getLociNames(g)[locus]
  file <- paste0(label, ".arp")
  
  # Header
  write("[Profile]", file = file)
  write(paste0("Title=\"", getDescription(g), "\""), file = file, append = TRUE)
  write(paste0("NbSamples=", getNumStrata(g)), file = file, append = TRUE)
  data.type <- if(getPloidy(g) > 1) {
    "MICROSAT"
  } else if(is.null(getSequences(g)[[locus]])) {
    "FREQUENCY"
  } else {
    "DNA"
  }
  write(paste0("DataType=", data.type), file = file, append = TRUE)
  g.data <- paste0("GenotypicData=", ifelse(getPloidy(g) == 1, 0, 1))
  write(g.data, file = file, append = TRUE)
  write("GameticPhase=0", file = file, append = TRUE)
  write("MissingData='?'", file = file, append = TRUE)
  write("LocusSeparator=WHITESPACE", file = file, append = TRUE)
  
  # Data
  write("[Data]", file = file, append = TRUE)
  write("[[Samples]]", file = file, append = TRUE)  
  for(st in strataSplit(g)) {
    write(
      paste0("SampleName=\"", getStrata(st)[1], "\""), 
      file = file, 
      append = TRUE
    )
    write(paste0("SampleSize=", getNumInd(st)), file = file, append = TRUE)
    write("SampleData={", file = file, append = TRUE)
    if(getPloidy(g) == 1) { # Sequences
      hap.df <- as.data.frame(alleleFreqs(st)[[locus]])
      colnames(hap.df)[1] <- locus
      dna <- getSequences(st)[[locus]]
      if(!is.null(dna)) {
        dna <- dna %>% 
          as.matrix() %>% 
          as.character() %>% 
          toupper()
        dna <- apply(dna, 1, paste, collapse = "")[hap.df[[locus]]]
      }
      hap.df <- cbind(hap.df, dna)
      write(
        t(hap.df), 
        ncolumns = ncol(hap.df), 
        sep = "\t", 
        file = file, 
        append = TRUE
      )
    } else { # Microsats
      loc.mat <- .stackedAlleles(st, na.val = "?") %>% 
        dplyr::mutate(
          id = ifelse(
            .data$allele == 1, 
            .data$id, 
            sapply(nchar(.data$id), function(i) {
              paste(rep(" ", i), collapse = "")
            })
          ),
          allele = ifelse(.data$allele == 1, "1", " ")
        ) %>% 
        dplyr::select(.data$id, .data$allele, dplyr::everything()) %>% 
        as.matrix()
      
      write(
        t(loc.mat), 
        ncolumns = ncol(loc.mat), 
        sep = "\t", 
        file = file, 
        append = TRUE
      )
    }    
    write("}", file = file, append = TRUE)
  }
  
  # Structure
  write("[[Structure]]", file = file, append = TRUE)
  st.name <- paste(
    "A group of", getNumStrata(g), "populations analyzed for", data.type
  )
  st.name <- paste0("StructureName=\"", st.name, "\"")
  write(st.name, file = file, append = TRUE)
  write("NbGroups=1", file = file, append = TRUE)
  write("Group= {", file = file, append = TRUE)
  for(st in getStrataNames(g)) {
    write(paste0("\"", st, "\""), file = file, append = TRUE)
  }
  write("}", file = file, append = TRUE)
  
  invisible(NULL)
}