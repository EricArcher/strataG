#' @title Hardy-Weinberg Equilibrium
#' @description Calculate Hardy-Weinberg equilibrium p-values.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param use.genepop logical. Use GENEPOP to calculate HWE p-values? 
#'   If \code{FALSE} then \code{\link[pegas]{hw.test}} is used.
#' @param which,enumeration,dememorization,batches parameters for GENEPOP MCMC
#'   HWE procedure as defined in \code{\link[genepop]{test_HW}}.
#' @param num.rep the number of replicates for the Monte Carlo procedure for
#'   \code{\link[pegas]{hw.test}} or number of iterations for
#'   \code{\link[genepop]{test_HW}}.
#' @param delete.files logical. Delete GENEPOP files when done?
#' @param label character string to use to label GENEPOP files.
#'   
#' @return a vector of p-values for each locus.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \code{\link{genepop}}, \code{\link[pegas]{hw.test}}
#' 
#' @examples
#' data(msats.g)
#' hweTest(msats.g)
#' 
#' @aliases hwe HWE
#' @export
#' 
hweTest <- function(
  g, 
  use.genepop = FALSE,
  which = c("Proba", "excess", "deficit"),
  enumeration = FALSE,
  dememorization = 10000,
  batches = 20, 
  num.rep = 5000, 
  delete.files = TRUE, 
  label = NULL
) {
  
  if(use.genepop) {
    g <- stratify(g)
    in.file <- genepopWrite(g, label)
    
    # Run Genepop  
    output <- genepop::test_HW(
      in.file$fname, 
      which = match.arg(which),
      enumeration = enumeration,
      dememorization = dememorization,
      batches = batches,
      iterations = num.rep,
      verbose = FALSE
    )
    result <- scan(output, what = "character", quiet = TRUE)

    # Find HWE p-value
    result[result == "-"] <- NA
    result[result == "No"] <- NA
    i <- which(result %in% names(in.file$locus.names))
    hwe.p <- as.numeric(result[i + 1])
    names(hwe.p) <- in.file$locus.names[result[i]]
    
    if(delete.files) {
      files <- c(output, in.file$fname, "fichier.in", "cmdline.txt")
      for(f in files) if(file.exists(f)) file.remove(f)
    }
    
    hwe.p
  } else {
    hwe <- pegas::hw.test(gtypes2genind(g), B = num.rep)
    hwe[, ncol(hwe)]
  }
}