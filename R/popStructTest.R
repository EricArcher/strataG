#' @name popStructTest 
#' @title Population Differentiation Tests
#' @description Conduct overall and/or pairwise tests of 
#'   population differentiation.
#' 
#' @param g a \code{\link{gtypes}} object.
#' @param nrep number specifying number of permutation replicates to use for 
#'   permutation test.
#' @param stats a character vector or list of functions specifying which
#'   anlayses to conduct. If characters, then valid possible choices are:
#'   "phist", "fst", "fst.prime", "fis", "gst", "gst.prime", "gst.dbl.prime",
#'   "d", or "chi2", or "all". If a list, then functions must be a valid
#'   population structure function (see \code{\link{popStructStat}}) taking a
#'   \linkS4class{gtypes} object and returning a named statistic estimate.
#' @param type character determining type of test to conduct. Can be "overall", 
#'   "pairwise", or "both". If "pairwise" or "both" are chosen and there are 
#'   only two strata, then only an overall test will be conducted.
#' @param max.cores The maximum number of cores to use to distribute separate 
#'   statistics over. The number of cores to use to distribute separate
#'   statistics over. If set to \code{NULL}, the value will be what is 
#'   reported by \code{\link[parallel]{detectCores} - 1}. 
#'   If \code{detectCores} reports \code{NA}, 
#'   \code{max.cores} will be set to 1.
#' @param keep.null logical. Keep the null distribution from the 
#'   permutation test?
#' @param quietly logical. Print progress and results?
#' @param write.output logical. Write a .csv file with results?
#' @param ... other parameters to be passed to population 
#'   differentiation functions.
#' 
#' @return \describe{
#'  \item{overall}{a list containing: \describe{
#'    \item{\code{strata.freq}}{a vector of the sample sizes for each stratum}
#'    \item{\code{result}}{a matrix with the statistic estimate and p-value 
#'      for each statistic}
#'    \item{\code{null.dist}}{a matrix with the null distributions for 
#'      each statistic}
#'  }}
#'  \item{pairwise}{a list containing: \describe{
#'    \item{\code{result}}{a data.frame with the result of each pairwise 
#'      comparison on each row}
#'    \item{\code{pair.mat}}{a list with a pairwise matrix for each statistic.
#'      Values in lower left are the statistic estimate, and upper right are
#'      p-values}
#'    \item{\code{null.dist}}{a matrix with the null distributions for 
#'      each statistic}
#'  }} 
#' }
#' 
#' @note On multi-core systems, runs of separate statistics are automatically 
#'   distributed over as many cores as available (minus one). This can be 
#'   controlled by the \code{max.cores} argument if less core usage is 
#'   desired. 
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples 
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#' 
#' # Just an overall Chi-squared test
#' ovl <- overallTest(msats.g, stats = "chi2", nrep = 100)
#' ovl
#' 
#' #' Just a pairwise test for Gst
#' pws <- pairwiseTest(msats.g, stats = "gst", nrep = 100)
#' pws
#' 
#' \dontrun{
#' #' Both overall and pairwise tests for Fst and F'st
#' full <- popStructTest(msats.g, stats = c("fst", "fst.prime"))
#' print(full$overall)
#' print(full$pairwise)
#' }
#' 
#' @export
#' 
popStructTest <- function(g, nrep = 1000, stats = "all", 
                          type = c("both", "overall", "pairwise"),
                          keep.null = FALSE, quietly = FALSE, 
                          max.cores = 1, write.output = FALSE, ...) {
  # check arguments
  type <- match.arg(type)
    
  # conduct overall test
  overall <- NULL
  if(type %in% c("both", "overall")) {
    overall <- overallTest(
      g = g, nrep = nrep, stats = stats, keep.null = keep.null, 
      quietly = quietly, max.cores = max.cores, ...
    )    
  }
  
  # conduct pairwise test
  pairwise <- NULL
  if(type %in% c("both", "pairwise") & getNumStrata(g) > 2) {
    pairwise <- pairwiseTest(
      g = g, nrep = nrep, stats = stats, keep.null = keep.null, 
      quietly = quietly, max.cores = max.cores, ...
    )
  }

  if(write.output) {
    if(!is.null(overall)) {
      out.file <- gsub("[[:punct:]]", ".", 
       paste(getDescription(g), "permutation test results.csv")
      )
      utils::write.csv(overall$result, out.file)
    }
    if(!is.null(pairwise)) {
      for(stat in names(pairwise$pair.mat)) {
        out.file <- gsub("[[:punct:]]", ".", 
          paste(getDescription(g), stat, "pairwise matrix.csv")
        )
        utils::write.csv(pairwise$pair.mat[[stat]], out.file)
      }
    }    
  }
    
  invisible(list(overall = overall, pairwise = pairwise))
}


#' @keywords internal
#' 
.checkStats <- function(stats, ploidy) {
  avail.stats <- if(ploidy == 1) {
    c("chi2", "fst", "phist")
  } else {
    c(
      "chi2", "d", 
      "fis", "fst", "fst.prime", 
      "gst", "gst.prime", "gst.dbl.prime"
    )
  }
  stats <- tolower(stats)
  stats <- if("all" %in% stats) {
    avail.stats 
  } else {
    missing <- setdiff(stats, avail.stats)
    if(length(missing) > 0) {
      missing <- paste(missing, collapse = ", ")
      stop(paste("the following stats are not available:", missing))
    }
    unique(stats)
  }
  if(length(stats) == 0) stop("no stats specified.")
  stats
}


#' @noRd
#' 
.runStatFunc <- function(stat.name, input) {
  stat.func <- switch(
    tolower(stat.name),
    chi2 = .statChi2,
    d = .statJostD,
    fis = .statFis,
    fst = .statFst,
    fst.prime = .statFstPrime,
    gst = .statGst,
    gst.prime = .statGstPrime,
    gst.dbl.prime = .statGstDblPrime,
    phist = .statPhist,
    ... = NULL
  )
  if(is.null(stat.func)) {
    .formatResult(stat.name, NULL, input$keep.null)
  } else {
    stat.func(input)
  }
}


#' @rdname popStructTest
#' @export
#' 
overallTest <- function(g, nrep = 1000, stats = "all", keep.null = FALSE, 
                        quietly = FALSE, max.cores = 1, ...) {
  # check requested stats
  stats <- .checkStats(stats, getPloidy(g))
  
  # check replicates
  if(is.null(nrep)) nrep <- 0
  if(!is.numeric(nrep) & length(nrep) != 1) {
    stop("'nrep' must be a single-element numeric vector")
  }
  if(nrep == 0) keep.null <- FALSE
  
  # format gtypes input
  input <- .formatCinput(g, nrep, keep.null, hap.dist = "phist" %in% stats, ...)
  # collect strata frequencies to named vector 
  strata.freq <- .strataFreq(g)
  desc <- getDescription(g)
  rm(g)
  
  # run population structure tests
  if(!quietly) cat(
    cat("\n<<<", desc, ">>>\n"),
    format(Sys.time()), ": Overall test :", nrep, "permutations\n"
  )
  cl <- swfscMisc::setupClusters(length(stats), max.cores)
  result <- tryCatch({
    if(is.null(cl)) {
      lapply(stats, .runStatFunc, input = input)
    } else {
      parallel::clusterEvalQ(cl, require(strataG))
      parallel::clusterExport(cl, "input", environment())
      parallel::parLapplyLB(cl, stats, .runStatFunc, input = input)
    }
  }, finally = if(!is.null(cl)) parallel::stopCluster(cl) else NULL)
  
  # create matrix of estimates and p-values
  result.mat <- t(sapply(result, function(x) x$result))
  rownames(result.mat) <- sapply(result, function(x) x$stat.name)
  
  # collect null distribution
  null.dist <- if(keep.null) {
    nd <- sapply(result, function(x) x$null.dist)
    colnames(nd) <- rownames(result.mat)
    nd
  } else NULL
  
  if(!quietly) {
    cat("\n")
    print(cbind(N = strata.freq))
    cat("\nPopulation structure results:\n")
    print(result.mat)
    cat("\n")
  }
  
  invisible(
    list(strata.freq = strata.freq, result = result.mat, null.dist = null.dist)
  )
}


#' @rdname popStructTest
#' @export
#' 
pairwiseTest <- function(g, nrep = 1000, stats = "all", keep.null = FALSE, 
                         quietly = FALSE, max.cores = 1, ...) { 
  
  if(getNumStrata(g) == 1) stop("'g' must have more than one stratum defined.")
  
  if(!quietly) cat(
    cat("\n<<<", getDescription(g), ">>>\n"),
    format(Sys.time()), ": Pairwise tests :", nrep, "permutations\n"
  )
  
  # create strata pairs
  strata.pairs <- .strataPairs(g)
  
  # run permutation test on all pairwise gtypes subsets
  pair.list <- vector("list", length = nrow(strata.pairs))
  for(i in 1:nrow(strata.pairs)) {
    pair <- unlist(strata.pairs[i, ])
    
    if(!quietly) {
      cat("  ", format(Sys.time()), ":", paste(pair, collapse = " v. "), "\n")
    }
    
    pair.list[[i]] <- overallTest(
      g = g[ , , pair], nrep = nrep, stats = stats, keep.null = keep.null, 
      quietly = TRUE, max.cores = max.cores, ...
    )
  }
  
  # compile results in 'pair.list' into a data.frame
  result <- do.call(rbind, lapply(pair.list, function(pair) {
    result.vec <- as.vector(t(pair$result))
    names(result.vec) <- paste(
      rep(rownames(pair$result), each = 2), c("", ".p.val"), sep = ""
    ) 
    result.vec <- rbind(result.vec)
    s1 <- names(pair$strata.freq)[1]
    s2 <- names(pair$strata.freq)[2]
    n1 <- pair$strata.freq[1]
    n2 <- pair$strata.freq[2]
    strata.1 <- paste(s1, " (", n1, ")", sep = "")
    strata.2 <- paste(s2, " (", n2, ")", sep = "")
    df <- data.frame(
      strata.1 = s1, 
      strata.2 = s2, 
      n.1 = n1, 
      n.2 = n2,
      stringsAsFactors = FALSE
    )
    df <- cbind(df, result.vec) 
    rownames(df) <- paste(strata.1, " v. ", strata.2, sep = "")
    df
  }))
  
  # create pairwise matrices - lower left is estimate, upper right is p-value 
  stat.cols <- seq(5, ncol(result), 2)
  strata <- getStrataNames(g)
  mat <- matrix(nrow = length(strata), ncol = length(strata), 
    dimnames = list(strata, strata)
  )
  pair.mat <- lapply(stat.cols, function(i) {
    for(j in 1:nrow(result)) {
      strata.1 <- as.character(result$strata.1[j])
      strata.2 <- as.character(result$strata.2[j])
      mat[strata.2, strata.1] <- result[j, i]
      mat[strata.1, strata.2] <- result[j, i + 1]
    }
    mat
  })
  names(pair.mat) <- colnames(result)[stat.cols]
  
  # compile null distributions into list of matrices
  null.dist <- if(keep.null) {
    null.mat <- lapply(pair.list, function(pair) pair$null.dist)
    names(null.mat) <- result$pair.label
    null.mat
  } else NULL
  
  if(!quietly) {
    cat("\nPopulation structure results:\n")
    print(result[, 5:ncol(result)])
    cat("\n")
  }
  
  invisible(list(result = result, pair.mat = pair.mat, null.dist = null.dist))
}