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
  result <- numericSNPmat(g)
  
  result[is.na(result)] <- "?"
  result <- lapply(1:nrow(result), function(i) result[i, ])
  strata <- gsub("[ _]", ".", strata(g))
  id <- gsub("[ _]", ".", indNames(g))
  names(result) <- paste(strata, id, sep = "_")
  
  write.nexus.data(result, file = file)
  
  snapp.file <- scan(file, what = "character", sep = "\n", quiet = TRUE)
  bgn <- grep("BEGIN", snapp.file)
  snapp.file[bgn] <- "BEGIN CHARACTERS;"
  fmt <- grep("FORMAT", snapp.file)
  snapp.file[fmt] <- "  FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS=\"012\" LABELS=LEFT TRANSPOSE=NO INTERLEAVE=NO;"
  write(snapp.file, file = file)
  
  invisible(do.call(rbind, result))
}