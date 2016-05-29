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
#' @aliases arlequin
#' @importFrom reshape2 dcast
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