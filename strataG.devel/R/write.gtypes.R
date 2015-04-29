#' @title Write gtypes
#' @description Write a \code{\link{gtypes}} object to file(s).
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param label label for filename(s). Default is the gtypes description 
#'   if present.
#' @param folder folder where file(s) should be written to. If NULL, files 
#'   are written to current working directory.
#' @param as.frequency logical indicating if haploid data should be output 
#'   as frequency tables.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
write.gtypes <- function(g, label = NULL, folder = NULL, as.frequency = FALSE) {
  desc <- description(g)
  label <- if(!is.null(label)) {
    label 
  } else if(!is.null(desc)) {
    desc 
  } else "strataG.gtypes"
  fname <- paste(label, ".csv", sep = "")
  if(!is.null(folder)) fname <- file.path(folder, fname)
  csv.file <- if(ploidy(g) == 1 & as.frequency) {
    x <- as.frequency(g) 
    x <- data.frame(haplotype = rownames(x), cbind(x))
    rownames(x) <- NULL
    x
  } else as.matrix(g)
  csv.file <- cbind(id = rownames(csv.file), strata = strata(g), csv.file)
  write.csv(csv.file, file = fname, row.names = FALSE)
  dna <- sequences(g)
  if(!is.null(dna)) {
    fname <- paste(label, ".fasta", sep = "")
    if(!is.null(folder)) fname <- file.path(folder, fname)
    write.fasta(dna, fname)
  }
}