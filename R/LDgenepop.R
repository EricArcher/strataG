#' @title Linkage Disequilibrium
#' @description Calculate linkage disequilibrium p-values using GENEPOP.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param dememorization,batches,iterations parameters for GENEPOP MCMC
#'   LD procedure as defined in \code{\link[genepop]{test_LD}}.
#' @param delete.files logical. Delete GENEPOP input and output files when done?
#' @param label character string to use to label GENEPOP input and output files.
#' 
#' @return data.frame of disequilibrium estimates between pairs of loci
#' @details Samples are pooled into a single population before testing;
#'   existing strata are not tested separately. At least two loci are required.
#'   Pairs reported by GENEPOP as having no contingency table receive
#'   \code{NA} for the p-value, standard error, and number of switches.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{genepop}
#' 
#' @examples \dontrun{
#' data(msats.g)
#' msats.ld <- LDgenepop(msats.g)
#' head(msats.ld)
#' }
#' 
#' @export
#' 
LDgenepop <- function(
  g, 
  dememorization = 10000,
  batches = 100, 
  iterations = 5000, 
  delete.files = TRUE, 
  label = NULL
) {
    
  if(getNumLoci(g) < 2L) {
    stop("At least two loci are required for linkage disequilibrium testing.",
         call. = FALSE)
  }
  # Run Genepop
  g <- stratify(g)
  in.file <- genepopWrite(g, label)
  
  output <- genepop::test_LD(
    in.file$fname,
    dememorization = dememorization,
    batches = batches,
    iterations = iterations,
    verbose = FALSE
  )
  # Keep three result fields for untestable pairs without coercion warnings.
  lines <- readLines(output, warn = FALSE)
  lines <- gsub("No[[:space:]]+contingency[[:space:]]+table",
                "NA NA NA", lines)
  result <- scan(text = paste(lines, collapse = "\n"),
                 what = "character", quiet = TRUE)
  
  # Create empty matrix
  locus.names <- in.file$locus.names
  numrows <- ((length(locus.names) ^ 2) - length(locus.names)) / 2
  result.mat <- matrix(
    as.character(NA), 
    nrow = numrows, 
    ncol = 5,
    dimnames = list(
      1:numrows,
      c("Locus.1", "Locus.2", "p.value", "std.err", "switches")
    )
  )
  
  # Find starting points
  loc <- grep("Switches", result, value = F) + 7
  first.col <- grep(names(locus.names)[1], result)[1]
  num.skip <- first.col - loc
  row.mask <- c(rep(F, num.skip), rep(T, 5))
  
  # Read matrix
  for(r in 1:numrows) {
    result.mat[r, ] <- result[loc:(loc + num.skip + 5)][row.mask]
    loc <- loc + num.skip + 5
  }
  
  # Convert to data.frame and format columns
  result.mat[result.mat == "NA"] <- NA_character_
  result.df <- tibble::as_tibble(result.mat) |> 
    dplyr::mutate(
      Locus.1 = locus.names[.data$Locus.1],
      Locus.2 = locus.names[.data$Locus.2],
      p.value = as.numeric(.data$p.value),
      std.err = as.numeric(.data$std.err),
      switches = as.integer(.data$switches)
    ) |> 
    as.data.frame()
  
  if(delete.files) {
    files <- c(output, in.file$fname, "fichier.in", "cmdline.txt")
    for(f in files) if(file.exists(f)) file.remove(f)
  }
  
  result.df
}
