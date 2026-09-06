#' @title Write NEXUS File for SNAPP
#' @description Write NEXUS File for SNAPP
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param file the filename the NEXUS file to output.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
write.nexus.snapp <- function(g, file = "snapp.data.nex") {
  result <- as.data.frame(g, coded.snps = TRUE)
  strata <- gsub("[ _]", ".", result$stratum)
  id <- gsub("[ _]", ".", result$id)
  result$id <- result$stratum <- NULL
  
  result <- lapply(seq_len(nrow(result)), function(i) {
    unlist(result[i, , drop = FALSE], use.names = TRUE)
  })
  names(result) <- paste(strata, id, sep = "_")
  
  file.data <- lapply(result, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- "?"
    x
  })
  ape::write.nexus.data(file.data, file = file, format = "standard",
                        interleaved = FALSE)
  
  snapp.file <- scan(file, what = "character", sep = "\n", quiet = TRUE)
  bgn <- grep("BEGIN", snapp.file)
  snapp.file[bgn] <- "BEGIN CHARACTERS;"
  fmt <- grep("FORMAT", snapp.file)
  snapp.file[fmt] <- "  FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS=\"012\" LABELS=LEFT TRANSPOSE=NO INTERLEAVE=NO;"
  write(snapp.file, file = file)
  
  invisible(do.call(rbind, result))
}
