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
#' writeGtypes(msats.g, one.col = TRUE)
#' 
#' # Write control region data as frequency tables
#' data(dloop.g)
#' writeGtypes(dloop.g, as.frequency = TRUE)
#' }
#' 
#' @export
#' 
writeGtypes <- function(g, label = NULL, folder = NULL, by.strata = TRUE,
                        as.frequency = FALSE, freq.type = c("freq", "prop"),
                        ...) {
  label <- .getFileLabel(g, label)
  
  g.mats <- if(getPloidy(g) == 1 & as.frequency) {
    alleleFreqs(g, by.strata = by.strata, type = freq.type) %>% 
      sapply(function(x) {
        as.data.frame.matrix(x) %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column("id") %>% 
          dplyr::select(.data$id, dplyr::everything())
      }, simplify = FALSE) %>% 
      stats::setNames(paste(label, names(.data)))
  } else {
    stats::setNames(list(as.matrix(g, ...)), label)
  }
  
  out.files <- NULL
  for(f in names(g.mats)) {
    fname <- paste0(f, ".csv")
    if(!is.null(folder)) fname <- file.path(folder, fname)
    out.files <- c(out.files, f)
    utils::write.csv(g.mats[[f]], file = fname, row.names = FALSE)
  }
  
  if(!is.null(getSequences(g))) {
    for(x in getLociNames(g)) {
      fname <- paste(label, x, "fasta", sep = ".")
      if(!is.null(folder)) fname <- file.path(folder, fname)
      ape::write.dna(
        getSequences(g)[[x]], 
        file = fname, 
        format = "fasta", 
        nbcol = -1, 
        colsep = "", 
        indent = 0, 
        blocksep = 0
      )
      out.files <- c(out.files, fname)
    }
  }
  
  invisible(out.files)
}