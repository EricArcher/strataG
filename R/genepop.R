#' @name genepop
#' @title Run GENEPOP
#' @description Format output files and run GENEPOP. Filenames used are returned 
#'   so that output files can be viewed or read and parsed into R.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param output.ext character string to use as extension for output files.
#' @param show.output logical. Show GENEPOP output on console?
#' @param label character string to use to label GENEPOP input and output files.
#' @param dem integer giving the number of MCMC dememorisation or burnin steps.
#' @param batches integer giving number of MCMC batches.
#' @param iter integer giving number of MCMC iterations.
#' @param other.settings character string of optional GENEPOP command line 
#'   arguments.
#' @param input.fname character string to use for input file name.
#' @param exec name of Genepop executable
#' 
#' @note GENEPOP is not included with \code{strataG} and must be downloaded 
#'   separately. Additionally, it must be installed such that it can be run from 
#'   the command line in the current working directory. See the vignette 
#'   for \code{external.programs} for installation instructions.
#' 
#' @return \describe{
#'   \item{\code{genepop}}{a list with a vector of the locus names and a vector 
#'     of the input and output filenames.}
#'   \item{\code{genepopWrite}}{a list with the filename written and a vector
#'     mapping locus names in the file to the original locus names.}
#' }
#' 
#' @references GENEPOP 4.3 (08 July 2014; Rousset, 2008)\cr
#'   \url{http://kimura.univ-montp2.fr/~rousset/Genepop.htm}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @seealso \link{hweTest}, \link{LDgenepop}
#' 
#' @examples \dontrun{
#' # Estimate Nm for the microsatellite data
#' data(msats.g)
#' # Run Genepop for Option 4
#' results <- genepop(msats.g, output.ext = ".PRI", other.settings = "MenuOptions=4")
#' # Locus name mapping and files
#' results
#' # Show contents of output file
#' file.show(results$files["output.fname"])
#' }
#' 
#' @export
#' 
genepop <- function(g, output.ext = "", show.output = F, label = "genepop.run",
                    dem = 10000, batches = 100, iter = 5000, 
                    other.settings = "", input.fname = "loc_data.txt",
                    exec = "Genepop") {
  
  locus.names <- genepopWrite(g, label)
  
  # Write settings file
  settings.fname <- "settings.txt"
  write(c(
    paste("InputFile=", input.fname, sep = ""),
    "Mode=Batch",
    paste("Dememorisation=", as.integer(ifelse(dem < 100, 100, dem)), sep = ""),
    paste("BatchNumber=", as.integer(ifelse(batches < 10, 10, batches)), sep = ""),
    paste("BatchLength=", as.integer(ifelse(iter < 400, 400, iter)), sep = ""),
    other.settings
  ), file = settings.fname)
  
  # Run Genepop
  args <- list(
    command = paste0(exec, " settingsFile=", settings.fname), 
    intern = F, 
    ignore.stdout = !show.output, 
    wait = T, 
    ignore.stderr = T
  )
  if(.Platform$OS.type == "windows") {
    args <- c(
      args, 
      list(show.output.on.console = F, minimized = F, invisible = T)
    )
  }
  err.code <- do.call(system, args)

  if(err.code == 127) {
    stop("You do not have Genepop installed.")
  } else if(err.code != 0) {
    stop(paste("Error code", err.code, "returned from Genepop.")) 
  } 
  if(show.output) cat("\n")
  files <- c(settings.fname = settings.fname, input.fname = input.fname, 
    output.fname = paste(input.fname, output.ext, sep = ""), 
    cmdline = "cmdline.txt", fichier = "fichier.in"
  )
  invisible(list(locus.names = locus.names, files = files))
}


#' @rdname genepop
#' @export
#' 
genepopWrite <- function(g, label = NULL) {
  if(getPloidy(g) != 2) stop("'g' must be a diploid object")

  # convert alleles to equal spaced, zero-padded numbers
  .convertAlleles <- function(x) {
    x <- as.numeric(factor(x))
    x[is.na(x)] <- 0
    max.width <- formatC(x, width = max(2, nchar(x)), flag = "0")
  }
  
  # convert alleles to genotypes in GENEPOP format, substitute , for _ in
  # stratum and id labels
  genepop.fmt <- g@data %>% 
    dplyr::group_by(.data$locus) %>% 
    dplyr::mutate(allele = .convertAlleles(.data$allele)) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(.data$stratum, .data$id, .data$locus) %>% 
    dplyr::summarize(genotype = paste(.data$allele, collapse = "")) %>% 
    dplyr::ungroup() %>% 
    tidyr::spread(.data$locus, .data$genotype) %>% 
    dplyr::mutate(
      stratum = gsub(",", "_", .data$stratum),
      id = gsub(",", "_", .data$id)
    ) %>% 
    dplyr::ungroup()

  # create new short locus names to replace current ones
  locus.names <- stats::setNames(
    colnames(genepop.fmt)[-(1:2)],
    paste("LOC", 1:(ncol(genepop.fmt) - 2), sep = "")
  )
  
  # write header and locus names
  fname <- paste0(.getFileLabel(g, label), "_loc_data.txt")
  fname <- gsub(" ", "_", fname)
  write(getDescription(g), file = fname)
  write(
    paste(names(locus.names), collapse = ", "),
    file = fname, 
    append = TRUE
  )
  
  # write population genotypes
  genepop.fmt <- split(genepop.fmt, genepop.fmt$stratum)
  for(pop in genepop.fmt) {
    write("POP", file = fname, append = TRUE)
    for(i in 1:nrow(pop)) {
      id <- paste(pop[i, 1:2], collapse = " ")
      genotypes <- paste(pop[i, -(1:2)], collapse = " ")
      write(
        paste(id, genotypes, sep = " , "),
        file = fname, 
        append = TRUE
      )
    }  
  }
  
  invisible(list(fname = fname, locus.names = locus.names))
}