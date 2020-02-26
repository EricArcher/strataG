#' @title Read and Write Arlequin Files
#' @description Read an Arlequin-formatted project input file (.arp). Convert 
#'   .arp data into \code{gtypes} object. Write an input file from a 
#'   \code{gtypes} object.
#' 
#' @param file filename of an arlequin project (.arp) file. See
#'   \code{Notes} for details on how .arp files are parsed.
#' @param arp a list of arlequin profile information and data as returned from 
#'   \code{arlequinRead}.
#' @param avoid.dups logical. Should sample identifiers be combined with 
#'   strata names to avoid duplicate identifiers between strata? If set to 
#'   \code{FALSE}, ids will be left unchanged, but an error will be thrown
#'   when the \code{gtypes} object is created if duplicated ids are found.
#' @param g a \linkS4class{gtypes} object.
#' @param locus numeric or character designation of which locus to write for 
#'   haploid data.
#'   
#' @note \code{arp2gtypes()} will not create a \code{gtypes} object for 
#'   Arlequin projects with relative frequency data (\code{DataType=FREQUENCY}
#'   and \code{FREQUENCY=REL}). If \code{DataType=DNA} and
#'   \code{GenotypicData=0}, sequences for each haplotype or individual are 
#'   assumed to be from a single locus.\cr 
#' 
#' @details 
#' \tabular{ll}{
#'   \code{arlequinRead} \tab parses a .arp file.\cr
#'   \code{arp2gtypes} \tab converts list from parsed .arp file to gtypes.\cr
#'   \code{arlequinWrite} \tab writes gtypes to .arp file.\cr
#' }
#' 
#' @return
#' \describe{
#'  \item{arlequinRead}{a list containing: \tabular{ll}{
#'    \code{file} \tab name and full path of .arp file that was read.\cr
#'    \code{profile.info} \tab list containing parameters in \code{[[Profile]]} 
#'      section of .arp file. All parameters are provided. Parameters unset in 
#'      .arp file are set to default values.\cr
#'    \code{data.info} \tab list containing data from \code{[[Data]]} section 
#'      of .arp file. Can contain elements for \code{haplotype.definition} 
#'      (a data.frame), \code{distance.matrix} (a matrix), \code{sample.data} 
#'      (a data.frame), or \code{genetic.structure} (a list).\cr
#'  }}
#'  \item{arp2gtypes}{a \linkS4class{gtypes} object.}
#'  \item{arlequinWrite}{the filename of the .arp file that was written.}
#' }
#' 
#' @references Excoffier, L.G. Laval, and S. Schneider (2005) 
#'   Arlequin ver. 3.0: An integrated software package for population genetics 
#'   data analysis. Evolutionary Bioinformatics Online 1:47-50.\cr
#'   Available at \url{http://cmpg.unibe.ch/software/arlequin3/}
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' # write test microsat data .arp file
#' f <- arlequinWrite(msats.g, tempfile())
#' 
#' # read .arp file and show structure
#' msats.arp <- arlequinRead(f)
#' print(str(msats.arp))
#' 
#' # convert parsed data to gtypes object
#' msats.arp.g <- arp2gtypes(msats.arp)
#' msats.arp.g
#' 
#' # compare to original
#' msats.g
#' 
#' @name arlequin
#' @aliases arlequin
#' 
NULL

#' @rdname arlequin
#' @export
#' 
arlequinRead <- function(file) {
  # read file
  arp <- scan(file, what = "character", sep = "\n", quiet = TRUE)
  arp <- stringi::stri_trim_both(arp)
  arp <- arp[arp != ""]
  
  # collect Profile information
  locus.separator <- .getValues("LocusSeparator", arp)
  gametic.phase <- as.numeric(.getValues("GameticPhase", arp, "[[:punct:]]+"))
  recessive.data <- as.numeric(.getValues("RecessiveData", arp, "[[:punct:]]+"))
  recessive.allele <- .getValues("RecessiveAllele", arp, "[[:punct:]]+")
  missing.data <- .getValues("MissingData", arp)
  frequency <- .getValues("Frequency", arp)
  frequency.threshold <- as.numeric(.getValues("FrequencyThreshold", arp))
  epsilon.value <- as.numeric(.getValues("EpsilonValue", arp))
  profile.info <- list(
    title = .getValues("Title", arp, "[[:punct:]]+"),
    nb.samples = as.numeric(.getValues("NbSamples", arp, "[[:punct:]]+")),
    data.type = .getValues("DataType", arp, "[[:punct:]]+"),
    genotypic.data = as.numeric(.getValues("GenotypicData", arp, "[[:punct:]]+")),
    locus.separator = ifelse(is.na(locus.separator), "WHITESPACE", locus.separator),
    gametic.phase = ifelse(is.na(gametic.phase), 1, gametic.phase),
    recessive.data = ifelse(is.na(recessive.data), 0, recessive.data),
    recessive.allele = ifelse(is.na(recessive.data), "null", recessive.allele),
    missing.data = ifelse(is.na(missing.data), "?", missing.data),
    frequency = ifelse(is.na(frequency), "ABS", frequency),
    frequency.threshold = ifelse(is.na(frequency.threshold), 1e-5, frequency.threshold),
    epsilon.value = ifelse(is.na(epsilon.value), 1e-7, epsilon.value)
  )
  
  # collect Data information
  data.info <- list()
  
  # Haplotype list
  hap.def.i <- which(arp == "[[HaplotypeDefinition]]")
  if(length(hap.def.i) == 1) {
    data.info$haplotype.definition$name <- .getValues("HaplListName", arp)
    hap.list <- .getExtern("HaplList", arp, folder = dirname(file))
    # parse rows by whitespace from .arp or EXTERN file
    hap.list <- if(is.null(hap.list)) {
      .extractDataBlock(.getLine("HaplList", arp), arp)
    } else {
      strsplit(hap.list, "[[:space:]]+")
    }
    data.info$haplotype.definition$list <- .parseHaplotypicData(
      hap.list, 
      "id", 
      profile.info$locus.separator,
      profile.info$data.type == "DNA"
    )
  }
  
  # Distance matrix
  dist.mat.i <- which(arp == "[[DistanceMatrix]]")
  if(length(dist.mat.i) == 1) {
    data.info$distance.matrix$name <- .getValues("MatrixName", arp)
    matrix.data <- .getExtern("MatrixData", arp, folder = dirname(file))
    matrix.data <- if(is.null(matrix.data)) {
      .extractDataBlock(.getLine("MatrixData", arp), arp)
    } else {
      matrix.data <- strsplit(matrix.data, "[[:space:]]+")
      matrix.data <- matrix.data[sapply(matrix.data, length) > 0]
      lapply(matrix.data, function(x) x[x != ""])
    }
    size <- length(matrix.data) - 1
    dist.mat <- matrix(NA, nrow = size, ncol = size)
    rownames(dist.mat) <- colnames(dist.mat) <- matrix.data[[1]]
    matrix.data <- matrix.data[-1]
    for(i in 1:size) {
      dist.row <- as.numeric(matrix.data[[i]])
      dist.mat[i, ] <- c(dist.row, rep(NA, size - length(dist.row)))
    }
    data.info$distance.matrix$data <- swfscMisc::copy.tri(dist.mat, "lower")
  }
  
  # Samples
  samples.i <- which(arp == "[[Samples]]")
  if(length(samples.i) == 1) {
    sample.name <- .getValues("SampleName", arp)
    sample.data <- lapply(.getLine("SampleData", arp), function(i) {
      sample.data.i <- .extractDataBlock(i, arp)
      # replace missing data with NA
      sample.data.i <- lapply(sample.data.i, function(x) {
        x[x == profile.info$missing.data] <- NA
        x
      })
      # parse genetic data
      if(profile.info$data.type == "FREQUENCY") {
        .parseFrequencyData(sample.data.i)
      } else if(profile.info$genotypic.data) {
        .parseGenotypicData(sample.data.i) 
      } else {
        .parseHaplotypicData(
          sample.data.i, 
          c("id", "freq"),
          profile.info$locus.separator,
          profile.info$data.type == "DNA"
        ) %>% 
          dplyr::mutate(freq = as.numeric(.data$freq))
      }
    })
    
    data.info$sample.data <- do.call(
      rbind, 
      mapply(
        function(name, size, data) {
          cbind(strata = name, data, stringsAsFactors = FALSE) %>% 
            dplyr::select(.data$id, dplyr::everything())
        },
        name = sample.name, data = sample.data,
        SIMPLIFY = FALSE, USE.NAMES = FALSE
      )
    )
  }
  
  # Genetic structure
  structure.i <- which(arp == "[[Structure]]")
  if(length(structure.i) == 1) {
    data.info$genetic.structure$name <- .getValues("StructureName", arp)
    data.info$genetic.structure$groups <- lapply(
      .getLine("Group", arp), 
      function(i) gsub("\"", "", unlist(.extractDataBlock(i, arp)))
    )
  }
  
  list(file = file, profile.info = profile.info, data.info = data.info)
}


#' @rdname arlequin
#' @export
#' 
arp2gtypes <- function(arp, avoid.dups = FALSE) {
  if(arp$profile.info$frequency == "REL") {
    stop("can't convert RELative FREQUENCY data to gtypes")
  }
  if(arp$profile.info$genotypic.data) {
    .diploid2gtype(arp, avoid.dups)
  } else {
    loc.sep <- switch(
      arp$profile.info$locus.separator,
      WHITESPACE = " ",
      TAB = " ",
      NONE = "",
      arp$profile.info$locus.separator
    )
    
    switch(
      arp$profile.info$data.type,
      FREQUENCY = .freq2gtype(arp, avoid.dups),
      DNA = .dna2gtype(arp, avoid.dups),
      RFLP = .standard2gtype(arp, avoid.dups, loc.sep),
      STANDARD = .standard2gtype(arp, avoid.dups, loc.sep)
    )
  }
}



#' @rdname arlequin
#' @export
#' 
arlequinWrite <- function(g, file = NULL, locus = 1) {
  folder <- if(is.null(file)) "." else dirname(file)
  if(!is.null(file)) file <- basename(file)
  file <- .getFileLabel(g, file)
  if(!grepl("(\\.arp)$", tolower(file))) file <- paste0(file, ".arp")
  file <- file.path(folder, file)
  if(is.numeric(locus)) locus <- getLociNames(g)[locus]
  
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
  loc.sep <- switch(
    data.type,
    MICROSAT = "WHITESPACE",
    FREQUENCY = "NONE",
    DNA = "NONE"
  )
  write(paste0("LocusSeparator=", loc.sep), file = file, append = TRUE)
  
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
        dplyr::select(-.data$stratum) %>% 
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
  
  invisible(file)
}



# Arlequin internal functions --------------------------------------------------

#' @noRd
#' 
.removeComments <- function(x) {
  gsub("[[:space:]]#[[:alnum:]|[:punct:]|[:space:]]+$", "", x)
}

#' @noRd
#'
# function to get value folowing "=" for a header title
.getValues <- function(title, arp, remove = NULL) {
  title <- paste0("^", title, "[[:space:]]*=[[:space:]]*")
  x <- grep(title, arp, ignore.case = TRUE, value = TRUE)
  if(length(x) == 0) return(NA) 
  x <- gsub(title, "", x)
  x <- .removeComments(x)
  x <- gsub("[']|\"", "", x)
  if(!is.null(remove)) x <- gsub(remove, "", x)
  stringi::stri_trim_both(x)
}

#' @noRd
#'
.getLine <- function(title, arp) {
  title <- paste0("^", title, "=")
  i <- grep(title, arp, ignore.case = TRUE)
  if(length(i) == 0) NULL else i
}

#' @noRd
#'
.getExtern <- function(title, arp, remove = NULL, folder = ".") {
  x <- .getValues(title, arp, remove)
  if(grepl("EXTERN", x)) {
    file <- file.path(folder, gsub("EXTERN[[:space:]]+", "", x))
    if(!file.exists(file)) stop("can't find EXTERN file: '", file, "'")
    scan(file, what = "character", sep = "\n", quiet = TRUE)
  } else NULL
}

#' @noRd
#'
.extractDataBlock <- function(start, arp) {
  ends <- grep("}", arp)
  end <- min(ends[ends > start])
  data.block <- .removeComments(arp[(start + 1):(end-1)])
  strsplit(data.block, "[[:space:]]+")
}

#' @noRd
#'
.parseFrequencyData <- function(sample.data) {
  data.df <- do.call(rbind, sample.data) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  if(ncol(data.df) != 2) {
    stop("'SampleData' must be 2 columns if 'DataType=FREQUENCY'")
  }
  data.df %>% 
    stats::setNames(c("id", "freq")) %>% 
    dplyr::mutate(freq = as.numeric(.data$freq))
}

#' @noRd
#'
.parseHaplotypicData <- function(data.block, id.cols, loc.sep, is.dna) {
  hap.data <- do.call(rbind, data.block)
  if(is.dna & ncol(hap.data) > length(id.cols)) {
    loc.data <- hap.data[, -(1:length(id.cols)), drop = FALSE]
    hap.data <- cbind(
      hap.data[, 1:length(id.cols)],
      apply(loc.data, 1, paste, collapse = "")
    )
  } else {
    has.1.locus.col <- ncol(hap.data) == length(id.cols) + 1
    if(has.1.locus.col & !loc.sep %in% c("WHITESPACE", "TAB")) {
      loc.sep <- ifelse(loc.sep == "NONE", "", loc.sep)
      loc.data <- do.call(rbind, strsplit(hap.data[, ncol(hap.data)], loc.sep))
      hap.data <- cbind(hap.data[, 1:length(id.cols)], loc.data)
    }
  }
  loc.names <- if(ncol(hap.data) > length(id.cols)) {
    paste0("locus_", 1:(ncol(hap.data) - length(id.cols)))
  } else NULL
  stats::setNames(
    as.data.frame(hap.data, stringsAsFactors = FALSE),
    c(id.cols, loc.names)
  )
}

#' @noRd
#'
.parseGenotypicData <- function(sample.data) {
  # sample.data is a list of character strings from a single SampleData block
  if(length(sample.data) %% 2 != 0) {
    stop("'SampleData' must have an even number of rows if 'GenotypicData=1'")
  }
  # create matrices of alternating rows
  allele.1 <- do.call(rbind, sample.data[c(T, F)])
  allele.2 <- do.call(rbind, sample.data[c(F, T)])
  if(ncol(allele.1) != ncol(allele.2) + 2) {
    stop("alternating rows in 'SampleData' must be two columns different") 
  }
  # remove id and frequency row from first allele matrix
  id.freq <- allele.1[, 1:2]
  colnames(id.freq) <- c("id", "freq")
  allele.1 <- allele.1[, -(1:2), drop = FALSE]
  # create data.frame combining matrices
  colnames(allele.1) <- colnames(allele.2) <- paste0("locus_", 1:ncol(allele.1))
  data.df <- rbind(
    cbind(id.freq, allele = 1, allele.1),
    cbind(id.freq, allele = 2, allele.2)
  ) %>% 
    as.data.frame(stringsAsFactors = FALSE) %>% 
    dplyr::mutate(
      freq = as.numeric(.data$freq),
      allele = as.numeric(.data$allele)
    )
  # sort by id (original order) and allele number
  data.df[order(match(data.df$id, id.freq[, "id"]), data.df$allele), ]
}

#' @noRd
#'
.expandByFreq <- function(x) {
  do.call(
    rbind,
    lapply(1:nrow(x), function(i) {
      freq <- x$freq[i]
      if(freq == 1) return(x[i, ])
      df.i <- x[rep(i, freq), ]
      df.i$id <- paste0(df.i$id, "_", 1:freq)
      df.i
    })
  )
}

#' @noRd
#'
.avoidDups <- function(x, avoid.dups) {
  if(avoid.dups) {
    if(anyDuplicated(x$id)) x$id <- paste0(x$strata, "_", x$id)
  }
  x
}

#' @noRd
#'
.diploid2gtype <- function(arp, avoid.dups) {    
  sample.data <- arp$data.info$sample.data %>% 
    tidyr::gather("locus", "value", dplyr::starts_with("locus_")) %>% 
    dplyr::mutate(locus = paste0(.data$locus, ".", .data$allele)) %>% 
    dplyr::select(-.data$allele) %>% 
    tidyr::spread(.data$locus, .data$value)
  if(any(sample.data$freq > 1)) sample.data <- .expandByFreq(sample.data)
  sample.data %>% 
    dplyr::select(-.data$freq) %>% 
    .avoidDups(avoid.dups) %>% 
    df2gtypes(ploidy = 2, description = arp$profile.info$title)
}

#' @noRd
#'
.freq2gtype <- function(arp, avoid.dups) {      
  g <- .expandByFreq(arp$data.info$sample.data) %>% 
    dplyr::select(-.data$freq) %>% 
    dplyr::mutate(hap = .data$id) %>% 
    .avoidDups(avoid.dups) %>% 
    df2gtypes(ploidy = 1, description = arp$profile.info$title)
  dist.mat <- arp$data.info$distance.matrix$data
  if(!is.null(dist.mat)) setOther(g, "dist.mat") <- dist.mat
  g
}

#' @noRd
#'
.dna2gtype <- function(arp, avoid.dups) {
  sample.data <- arp$data.info$sample.data
  haps <- arp$data.info$haplotype.definition$list
  if(is.null(haps)) haps <- sample.data[, c("id", "locus_1")]
  sample.data$locus_1 <- sample.data$id
  g <- sample.data %>% 
    .expandByFreq() %>% 
    dplyr::select(.data$id, .data$strata, .data$locus_1) %>%
    .avoidDups(avoid.dups) %>% 
    df2gtypes(
      ploidy = 1,
      sequences = ape::as.DNAbin(
        stats::setNames(strsplit(haps$locus_1, ""), haps$id)
      ),
      description = arp$profile.info$title
    )
  dist.mat <- arp$data.info$distance.matrix$data
  if(!is.null(dist.mat)) setOther(g, "dist.mat") <- dist.mat
  g
}

#' @noRd
#'
.standard2gtype <- function(arp, avoid.dups, loc.sep) {
  sample.data <- arp$data.info$sample.data
  haps <- arp$data.info$haplotype.definition$list
  if(!is.null(haps)) {
    haps <- data.frame(
      id = haps[, 1], 
      hap = apply(haps[, -1, drop = FALSE], 1, paste, collapse = loc.sep),
      stringsAsFactors = FALSE
    )
    sample.data <- dplyr::left_join(sample.data, haps, by = "id")
  } else if(nrow(sample.data) > 3) {
    sample.data <- cbind(
      sample.data[, 1:3],
      hap = apply(
        sample.data[, -(1:3), drop = FALSE], 
        1, 
        paste, 
        collapse = loc.sep
      ),
      stringsAsFactors = FALSE
    )
  } else {
    sample.data$hap <- sample.data$id
  }
  
  g <- sample.data %>% 
    .expandByFreq() %>% 
    dplyr::select(.data$id, .data$strata, .data$hap) %>% 
    .avoidDups(avoid.dups) %>% 
    df2gtypes(ploidy = 1, description = arp$profile.info$title)
  dist.mat <- arp$data.info$distance.matrix$data
  if(!is.null(dist.mat)) setOther(g, "dist.mat") <- dist.mat
  g
}


# Deprecated -------------------------------------------------------------------

#' @rdname arlequin
#' @export
#' 
read.arlequin <- function(file) {
  warning("'read.arlequin()' has been deprecated. use 'arlequinRead()' instead.")
  arlequinRead(file)
}

#' @rdname arlequin
#' @export
#' 
write.arlequin <- function(g, file = NULL, locus = 1) {
  warning("'write.arlequin()' has been deprecated. use 'arlequinWrite()' instead.")
  arlequinWrite(g, file, locus)
}