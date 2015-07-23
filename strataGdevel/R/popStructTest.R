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
#'   "phist", "fst", "fst.prime", "gst", "gst.prime", "gst.dbl.prime", "d", 
#'   or "chi2", or "all". If a list, then functions must be a valid population 
#'   structure function (see \code{\link{popStructStat}}) taking a 
#'   \linkS4class{gtypes} object and returning a named statistic estimate.
#' @param type character determining type of test to conduct. Can be "overall", 
#'   "pairwise", or "both". If "pairwise" or "both" are chosen and there are 
#'   only two strata, then only an overall test will be conducted.
#' @param keep.null logical. Keep the null distribution from the 
#'   permutation test?
#' @param quietly logical. Print progress and results?
#' @param num.cores number of CPU cores to use. Value is passed to 
#'   \code{\link[parallel]{mclapply}}.
#' @param write.output logical. Write a .csv file with results?
#' @param ... other parameters to be passed to population 
#'   differentiation functions.
#' 
#' @return
#' \describe{
#'  \item{overall}{a list containing:
#'    \tabular{ll}{
#'      \code{strata.freq} \tab a vector of the sample sizes for each stratum.\cr
#'      \code{result} \tab a matrix with the statistic estimate and p-value 
#'        for each statistic.\cr
#'      \code{null.dist} \tab a matrix with the null distributions for 
#'        each statistic.\cr
#'    }}
#'  \item{pairwise}{a list containing:
#'    \tabular{ll}{
#'      \code{result} \tab a data.frame with the result of each pairwise 
#'        comparison on a line.\cr
#'      \code{pair.mat} \tab a list with a pairwise matrix for each statistic. 
#'        Values in lower left are p-values, upper right is statistic estimate.\cr
#'      \code{null.dist} \tab a matrix with the null distributions for 
#'        each statistic.\cr
#'    }}
#' }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(dolph.msats)
#' data(dolph.strata)
#' msats.merge <- merge(dolph.strata[, c("ids", "fine")], dolph.msats, all.y = TRUE)
#' msats <- df2gtypes(msats.merge, ploidy = 2)
#' 
#' # Conduct an overall Chi-squared test
#' ovl <- overallTest(msats, nrep = 5, 
#'   stat.list = statList("chi2"), quietly = FALSE
#' )
#' 
#' # Conduct a pairwise test for Gst
#' pws <- pairwiseTest(msats, nrep = 5, 
#'   stat.list = list(statGst), quietly = FALSE
#' )
#' 
#' # Conduct both overall and pairwise tests for Fst and F'st
#' full <- popStructTest(msats, nrep = 5, stats = "fst")
#' print(full$overall)
#' print(full$pairwise)
#' 
#' @export
#' 
popStructTest <- function(g, nrep = 100, stats = "all", 
                          type = c("both", "overall", "pairwise"),
                          keep.null = FALSE, quietly = FALSE, num.cores = 1, 
                          write.output = FALSE, ...) {
  # check arguments
  type <- match.arg(type)
    
  # conduct overall test
  overall <- NULL
  if(type %in% c("both", "overall")) {
    overall <- overallTest(g = g, nrep = nrep, stats = stats, 
      keep.null = keep.null, num.cores = num.cores, quietly = quietly, ...
    )    
  }
  
  # conduct pairwise test
  pairwise <- NULL
  if(type %in% c("both", "pairwise") & nStrata(g) > 2) {
    pairwise <- pairwiseTest(g = g, nrep = nrep, stats = stats, 
      keep.null = keep.null, num.cores = num.cores, quietly = quietly, ...
    )
  }

  if(write.output) {
    if(!is.null(overall)) {
      out.file <- gsub("[[:punct:]]", ".", 
       paste(description(g), "permutation test results.csv")
      )
      write.csv(overall$result, out.file)
    }
    if(!is.null(pairwise)) {
      for(stat in names(pairwise$pair.mat)) {
        out.file <- gsub("[[:punct:]]", ".", 
          paste(description(g), stat, "pairwise matrix.csv")
        )
        write.csv(pairwise$pair.mat[[stat]], out.file)
      }
    }    
  }
    
  invisible(list(overall = overall, pairwise = pairwise))
}


#' @rdname popStructTest
#' @importFrom parallel makeForkCluster
#' @importFrom parallel parLapply 
#' @importFrom parallel stopCluster
#' @importFrom swfscMisc pVal
#' @export
#' 
overallTest <- function(g, nrep = 100, stats = "all", 
                        keep.null = FALSE, quietly = FALSE, 
                        num.cores = 1, ...) {  
  
  stat.list <- statList(stats)
  
  # check replicates
  if(!is.numeric(nrep) & length(nrep) != 1) {
    stop("'nrep' must be a single-element numeric vector")
  }
  if(nrep < 1) keep.null = FALSE
  
  # remove unstratified samples
  if(any(is.na(strata(g)))) g <- g[, , strataNames(g)]
  
  if(!quietly) cat(
    cat("\n<<<", description(g), ">>>\n"),
    format(Sys.time()), ": Overall test :", nrep, "permutations\n"
  )
  
  # calculate list of observed values for each population structure function
  result <- matrix(nrow = length(stat.list), ncol = 2)
  rownames(result) <- 1:nrow(result)
  colnames(result) <- c("estimate", "p.val")
  for(i in 1:nrow(result)) {
    x <- stat.list[[i]](g, ...)
    result[i, "estimate"] <- x
    rownames(result)[i] <- names(x)
  }
  
  # remove results and stats where estimate is NA
  to.remove <- which(apply(result, 1, function(x) is.na(x["estimate"])))
  if(length(to.remove) > 0) {
    result <- result[-to.remove, , drop = FALSE]
    stat.list <- stat.list[-to.remove]
  }    
  
  # conduct permutation test
  null.dist <- NULL
  if(nrep > 0 & length(stat.list) > 0) {
    st <- strata(g)  
    perm.func <- function(i, st, stat.list, g) {
      ran.strata <- sample(st)
      sapply(1:length(stat.list), function(j) {
        stat.list[[j]](g, strata = ran.strata, ...)
      })
    }
    
    if(num.cores > 1) {
      # setup clusters
      cl <- makeForkCluster(num.cores)
      tryCatch({
        # calculate matrix of null distributions
        null.dist <- do.call(
          rbind, 
          parLapply(cl, 1:nrep, perm.func, st = st, stat.list = stat.list, g = g)
        )
      }, finally = stopCluster(cl))
    } else {
      null.dist <- do.call(
        rbind, 
        lapply(1:nrep, perm.func, st = st, stat.list = stat.list, g = g)
      )
    }
    colnames(null.dist) <- rownames(result)
    
    # calculate vector of p-values
    for(x in rownames(result)) {
      result[x, "p.val"] <- pVal(result[x, "estimate"], null.dist[, x])
    }
  } 
  if(!keep.null) null.dist <- NULL
  
  # collect strata frequencies to named vector 
  strata.freq <- table(strata(g), useNA = "no")
  
  if(!quietly) {
    cat("\n")
    print(cbind(N = strata.freq))
    cat("\nPopulation structure results:\n")
    print(result)
    cat("\n")
  }
  
  invisible(
    list(strata.freq = strata.freq, result = result, null.dist = null.dist)
  )
}


#' @rdname popStructTest
#' @export
#' 
pairwiseTest <- function(g, nrep = 100, stats = "all", 
                         keep.null = FALSE, quietly = FALSE,
                         num.cores = 1, ...) { 
  
  if(!quietly) cat(
    cat("\n<<<", description(g), ">>>\n"),
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
    
    pair.list[[i]] <- overallTest(g = subset(g, strata = pair), nrep = nrep,
      stats = stats, keep.null = keep.null, quietly = TRUE, 
      num.cores = num.cores, ...
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
    df <- data.frame(pair.label = paste(strata.1, " v. ", strata.2, sep = ""), 
      strata.1 = s1, strata.2 = s2, n.1 = n1, n.2 = n2,
      stringsAsFactors = FALSE
    )
    cbind(df, result.vec)   
  }))
  rownames(result) <- NULL
  
  # create pairwise matrices - lower left is estimate, upper right is p-value 
  stat.cols <- seq(6, ncol(result), 2)
  strata <- sort(levels(strata(g)))
  mat <- matrix(nrow = length(strata), ncol = length(strata), 
    dimnames = list(strata, strata)
  )
  pair.mat <- lapply(stat.cols, function(i) {
    for(j in 1:nrow(result)) {
      strata.1 <- result$strata.1[j]
      strata.2 <- result$strata.2[j]
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
    print(result[, c(1, 6:ncol(result))])
    cat("\n")
  }
  
  invisible(list(result = result, pair.mat = pair.mat, null.dist = null.dist))
}


#' @rdname popStructTest
#' @export
#' 
statList <- function(stats = "all") {
  # check stats and return list of functions
  stat.list <- list(
    chi2 = statChi2,
    d = statJostD,
    fst = statFst,
    fst.prime = statFstPrime,
    gst = statGst,
    gst.prime = statGstPrime,
    gst.dbl.prime = statGstDblPrime,
    phist = statPhist
  )
  
  if(is.character(stats)) {
    stats <- tolower(stats)
    if("all" %in% stats) {
      stat.list 
    } else {
      missing <- !stats %in% names(stat.list)
      if(sum(missing) > 0) {
        missing <- paste(stats[missing], collapse = ", ")
        stop(paste("the following stats could not be found:", missing))
      }
      stat.list[stats]
    } 
  } else if(is.list(stats) & all(sapply(stats, is.function))) {
    stats
  } else {
    stop("'stats' is not a list of functions.")
  }
  
}