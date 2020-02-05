#' @name maverickRun
#' @title Run MavericK
#' @description Run MavericK clustering algorithm
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param params a list specifying parameters for MavericK. All parameters are 
#'   available and can be specified by partial matching. The function will automatically 
#'   specify parameters related to data formatting (data, headerRow_on, 
#'   missingData, ploidy, ploidyCol_on, popCol_on), so those will be ignored. 
#'   For a full list of available parameters and their definitions, see the 
#'   MavericK documentation distributed with the program.
#' @param label folder where input and output files will be written to.
#' @param data_fname file name of data input file.
#' @param param_fname file name of parameters file.
#' @param exec name of executable for MavericK.
#' 
#' @note MavericK is not included with \code{strataG} and must be downloaded 
#'   separately. It can be obtained from \url{http://www.bobverity.com/}. 
#'   Additionally, it must be installed such that it can be run from 
#'   the command line in the current working directory. See the vignette 
#'   for \code{external.programs} for OS-specific installation instructions.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @references 
#' Robert Verity and Richard Nichols. (2016) Estimating the number of 
#' subpopulations (K) in structured populations. Genetics \cr
#' Robert Verity and Richard Nichols. (2016) Documentation for MavericK 
#' software: Version 1.0
#' 
#' @export
#' 
maverickRun <- function(g, params = NULL, label = "MavericK_files",
                        data_fname = "data.txt", 
                        param_fname = "parameters.txt",
                        exec = "Maverick1.0.5") {
  
  if(getPloidy(g) < 2) stop("ploidy of 'g' must be 2 or greater")

  # set parameters
  p <- .maverickDefaultParams()
  if(!is.null(params)) {
    if(!is.list(params)) stop("'params' must be a list")
    i <- pmatch(names(params), names(p), nomatch = NULL)
    nomatch <- which(is.na(params))
    if(length(nomatch) > 0) {
      nomatch <- paste(names(params)[nomatch], collapse = ", ")
      stop(paste("Can't find the following parameters:", nomatch))
    }
    p[i] <- params
  }
  p$data <- data_fname
  p$headerRow_on <- TRUE
  p$missingData <- -9
  p$ploidy <- getPloidy(g)
  p$ploidyCol_on <- FALSE
  p$popCol_on <- getNumStrata(g) > 1
  
  # set directory
  wd <- getwd()
  unlink(label, force = TRUE)
  if(!dir.exists(label)) dir.create(label)
  setwd(label)
  
  # write parameters
  write(paste("Description:", getDescription(g)), file = param_fname)
  write(
    paste0("Created on: ", format(Sys.time()), "\n"), 
    file = param_fname, 
    append = TRUE
  )
  p <- p[order(names(p))]
  for(x in names(p)) {
    write(paste0(x, "\t", p[[x]]), file = param_fname, append = TRUE)
  }
  
  # write data
  df <- as.data.frame(g@data)
  colnames(df) <- make.names(colnames(df))
  df$ids <- make.names(df$ids)
  df$strata <- make.names(df$strata)
  for(i in 3:ncol(df)) df[[i]] <- as.numeric(df[[i]])
  df[is.na(df)] <- -9
  df <- df[order(df$strata, nchar(df$ids), df$ids), ]
  rownames(df) <- NULL
  utils::write.table(
    df, file = data_fname, quote = FALSE, sep = "\t", row.names = FALSE
  )
  
  # run MavericK
  cmd <- paste0(exec, " -parameters ", param_fname)
  err.code <- system(cmd)
  if(err.code == 127) {
    stop("You do not have MavericK properly installed.")
  } else if(!err.code == 0) {
    stop("Error running MavericK. Error code ", err.code, " returned.")
  }
  
  setwd(wd)
}

#' @keywords internal
#' 
.maverickDefaultParams <- function() {
  list(
    admix_on = FALSE,
    alpha = 1,
    alphaPropSD = 0.1,
    EMalgorithm_on = FALSE,
    EMiterations = 100,
    EMrepeats = 10,
    exhaustive_on = FALSE,
    fixAlpha_on = TRUE,
    fixLabels_on = TRUE,
    Kmax = 2,
    Kmin = 1,
    mainBurnin = 100,
    mainRepeats = 1,
    mainSamples = 1000,
    mainThinning = 1,
    outputComparisonStatistics = "outputComparisonStatistics.csv",
    outputComparisonStatistics_on = FALSE,
    outputEvanno = "outputEvanno.csv",
    outputEvanno_on = FALSE,
    outputEvidence = "outputEvidence.csv",
    outputEvidence_on = TRUE,
    outputEvidenceDetails = "outputEvidenceDetails.csv",
    outputEvidenceDetails_on = FALSE,
    outputEvidenceNormalised = "outputEvidenceNormalised.csv",
    outputEvidenceNormalised_on = TRUE,
    outputLikelihood = "outputLikelihood.csv",
    outputLikelihood_on = FALSE,
    outputLog = "outputLog.txt",
    outputLog_on = TRUE,
    outputMaxLike_admixFreqs = "outputMaxLike_admixFreqs.csv",
    outputMaxLike_admixFreqs_on = FALSE,
    outputMaxLike_alleleFreqs = "outputMaxLike_alleleFreqs.csv",
    outputMaxLike_alleleFreqs_on = FALSE,
    outputPosteriorGrouping = "outputPosteriorGrouping.csv",
    outputPosteriorGrouping_on = FALSE,
    outputQmatrix_gene = "outputQmatrix_gene.csv",
    outputQmatrix_gene_on = FALSE,
    outputQmatrix_ind = "outputQmatrix_ind.csv",
    outputQmatrix_ind_on = TRUE,
    outputQmatrix_pop = "outputQmatrix_pop.csv",
    outputQmatrix_pop.on = FALSE,
    outputQmatrix_structureFormat_on = FALSE,
    outputQmatrixError_gene = "outputQmatrixError_gene.csv",
    outputQmatrixError_gene_on = FALSE,
    outputQmatrixError_ind = "outputQmatrixError_ind.csv",
    outputQmatrixError_ind_on = FALSE,
    outputQmatrixError_pop = "outputQmatrixError_pop.csv",
    outputQmatrixError_pop_on = FALSE,
    suppressWarning1_on = FALSE,
    thermodynamic_on = TRUE,
    thermodynamicBurnin = 100,
    thermodynamicRungs = 21,
    thermodynamicSamples = 1000,
    thermodynamicThinning = 1
  )
}