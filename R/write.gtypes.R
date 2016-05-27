#' @title Write \code{gtypes}
#' @description Write a \linkS4class{gtypes} object to file(s).
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param label label for filename(s). Default is the gtypes description 
#'   if present.
#' @param folder folder where file(s) should be written to. If \code{NULL}, 
#'   files are written to current working directory.
#' @param as.frequency logical indicating if haploid data should be output 
#'   as frequency tables.
#' @param by.strata if \code{as.frequency == TRUE}, calculate frequencies by strata?
#' @param freq.type if \code{as.frequency == TRUE}, write absolute frequencies 
#'   (\code{"freq"}) or proportions (\code{"prop"}).
#' @param ... optional arguments controlling what information is included in the 
#'   genotype file and how it is formatted passed to \link[strataG]{as.matrix}.
#' 
#' @details Writes a comma-delimited (.csv) file of genotypes and if sequences 
#'   are present, a .fasta file for each locus. If haploid and \code{as.frequency} 
#'   is \code{TRUE}, then frequency tables for each locus are written to 
#'   separate files.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples \dontrun{
#' # Write microsatellites with one column per locus
#' data(msats.g)
#' write.gtypes(msats.g, one.col = TRUE)
#' 
#' # Write control region data as frequency tables
#' data(dloop.g)
#' write.gtypes(dloop.g, as.frequency = TRUE)
#' }
#' 
#' @importFrom utils write.csv
#' @export
#' 
write.gtypes <- function(g, label = NULL, folder = NULL, as.frequency = FALSE, 
                         by.strata = TRUE, freq.type = c("freq", "prop"), ...) {
  label <- .getFileLabel(g, label)
  
  g.mats <- if(ploidy(g) == 1 & as.frequency) {
    freq.list <- alleleFreqs(g, by.strata = by.strata) 
    freq.type <- match.arg(freq.type)
    x <- if(by.strata) {
      lapply(freq.list, function(x) x[, freq.type, ])
    } else {
      lapply(freq.list, function(x) x[, freq.type])
    }
    names(x) <- paste(label, names(x), sep = ".")
    x
  } else {
    x <- list(as.matrix(g, ...))
    names(x) <- label
    x
  }
  names(g.mats) <- paste(names(g.mats), ".csv", sep = "")
  if(!is.null(folder)) names(g.mats) <- file.path(folder, names(g.mats))
  for(f in names(g.mats)) write.csv(g.mats[[f]], file = f, row.names = FALSE)
  out.files <- names(g.mats)
  
  if(!is.null(sequences(g))) {
    for(x in locNames(g)) {
      fname <- paste(label, x, "fasta", sep = ".")
      if(!is.null(folder)) fname <- file.path(folder, fname)
      write.dna(sequences(g, x), file = fname, format = "fasta", nbcol = -1, 
                colsep = "", indent = 0, blocksep = 0)
      out.files <- c(out.files, fname)
    }
  }
  
  invisible(out.files)
}