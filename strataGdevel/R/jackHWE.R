#' @title Hardy-Weinberg Equlibrium Jackknife
#' @description Test influence of samples on Hardy-Weinberg equilibrium via 
#'   jackknife.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param exclude.num Number of samples to exclude at a time.
#' @param min.hwe.samples minimum samples needed to calculate HWE.
#' @param show.progress logical. Show progress of jackknife?
#' @param ... other arguments to be passed to GENEPOP.
#' @param jack.result result from run of \code{jackHWE}.
#' @param alpha critical value to determine if exclusion is "influential".
#' @param x result from a call to \code{jackInfluential}.
#' @param main main title for influential sample plots from 
#'   \code{plot.jack.influential}.
#' 
#' @details \tabular{ll}{
#'   \code{jack.hwe} \tab performs a HWE jackknife where all combinations 
#'     of \code{exclude.num} samples are left out and HWE is recalculated.\cr
#'   \code{jack.influential} \tab calculates odds.ratios between jackknife 
#'     HWE and observed HWE and identifies "influential" samples. Samples 
#'     are "influential" if the observed HWE p-value is < \code{alpha}, but 
#'     is > \code{alpha} when the samples are not present.\cr
#'   \code{plot.jack.influential} \tab creates a cumulative frequency plot 
#'     of all odds-ratios from \code{jack.influential}. A vertical dashed 
#'     line marks the smallest influential exclusion.\cr
#' }
#' 
#' @return \code{jack.hwe} returns a list with:
#' \item{obs}{a named vector of HWE p-values for each locus.}
#' \item{jack}{a \code{data.frame} of HWE p-values where each row is an 
#'   exclusion and columns are loci.}
#' \item{gtypes}{the original \code{gtypes} object.}\cr
#' \code{jack.influential} returns a list with:
#' \item{influential}{a \code{data.frame} of influential exclusions.}
#' \item{allele.freqs}{a \code{data.frame} listing the allele frequencies of 
#'   influential exclusions.}
#' \item{odds.ratio}{a \code{matrix} of odds ratios between exclusions (rows) 
#'   and loci (columns).}
#' 
#' @references Morin, P.A., R.G. LeDuc, F.I. Archer, K.K. Martien, 
#'   R. Huebinger, J.W. Bickham, and B.L. Taylor. 2009. Significant deviations 
#'   from Hardy-Weinberg equilibirum caused by low levels of microsatellite 
#'   genotyping errors. Molecular Ecology Resources 9:498-504.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @export
#' 
jackHWE <- function(g, exclude.num = 1, min.hwe.samples = 5, 
                    show.progress = TRUE, ...) {  
  
  if((nrow(g$genotypes) - exclude.num) < min.hwe.samples){
    stop(paste("'exclude.num' or 'min.HWE.samples' is too large to analyze this data"))
  }
  
  # Setup array of samples to exclude
  exclude.arr <- combn(indNames(g), exclude.num)
  
  # Calculate observed HWE and run jackknife 
  start.time <- Sys.time()
  nsteps <- ncol(exclude.arr) + 1 
  if(show.progress) {
    cat(format(start.time, "%Y-%m-%d %H:%M:%S"), "|", 1, "/", nsteps, "\n")
  }
  obs <- HWEgenepop(g, show.output = FALSE, ...)
  jack <- sapply(1:ncol(exclude.arr), function(i) {           
    if(show.progress) {
      now <- Sys.time()    
      avg.time <- difftime(now, start.time, units = "secs") / (i + 1)
      est.complete.time <- now + (avg.time * (nsteps - i))
      cat(format(now, "%Y-%m-%d %H:%M:%S"), "|", i + 1, "/", nsteps, "| ETA =", 
          format(est.complete.time, "%Y-%m-%d %H:%M:%S\n")
      )
    }
    to.keep <- setdiff(indNames(g), exclude.arr[, i])
    jack.gtypes <- g[to.keep, , ]
    HWEgenepop(jack.gtypes, show.output = FALSE, ...)
  })

  exclude.vec <- apply(exclude.arr, 2, function(x) paste(x, collapse = ", "))
  result <- list(
    obs = obs, 
    jack = data.frame(excluded = exclude.vec, t(jack), stringsAsFactors = FALSE), 
    gtypes = g
  )
  class(result) <- c("jack.result", class(result))
  result
}


#' @rdname jackHWE
#' @importFrom swfscMisc odds
#' @export
#' 
jackInfluential <- function(jack.result, alpha = 0.05) {  
  if(!inherits(jack.result, "jack.result")) {
    stop("'jack.result' is not from 'jack.hwe'")
  }
  
  obs <- jack.result$obs
  exclude.vec <- jack.result$jack$excluded
  jack <- t(as.matrix(jack.result$jack[, -1]))
  
  # Find influential samples and calculate odds and odds ratios
  obs.arr <- matrix(obs, length(obs), length(exclude.vec))
  which.infl <- which((obs.arr <= alpha) & (jack > alpha), arr.ind = TRUE)
  influential <- if(nrow(which.infl) > 0) {
    infl <- data.frame(excluded = exclude.vec[which.infl[, "col"]], 
                       stringsAsFactors = FALSE) 
    infl$locus <- names(obs)[which.infl[, "row"]]
    infl$obs.pval <- obs[which.infl[, "row"]]
    infl$jack.pval <- apply(which.infl, 1, function(i) jack[i[1], i[2]])
    infl$obs.odds <- odds(infl$obs.pval)
    infl$jack.odds <- odds(infl$jack.pval)
    infl$odds.ratio <- infl$jack.odds / infl$obs.odds
    infl <- infl[order(infl$odds.ratio, decreasing = TRUE), ] 
    rownames(infl) <- 1:nrow(infl) 
    infl  
  } else NULL
  
  # Calculate allele frequencies of influential samples 
  allele.freqs <- if(!is.null(influential)) {
    samples.loci <- unique(do.call(rbind, lapply(1:nrow(influential), function(i) {
      samples <- unlist(strsplit(influential$excluded[i], ", "))
      data.frame(id = samples, locus = influential$locus[i], 
                 stringsAsFactors = FALSE)
    })))
    formatted.freqs <- alleleFreqFormat(samples.loci, jack.result$gtypes)
    rownames(formatted.freqs) <- 1:nrow(formatted.freqs)
    colnames(formatted.freqs)[1] <- "id"
    formatted.freqs
  } else NULL
  
  # Calculate full odds ratio matrix
  odds.ratio <- t(odds(jack) / odds(obs.arr))
  rownames(odds.ratio) <- exclude.vec
  colnames(odds.ratio) <- names(obs)
  
  result <- list(influential = influential, allele.freqs = allele.freqs, 
                 odds.ratio = odds.ratio)
  class(result) <- c("jack.influential", class(result))
  result
}


#' @rdname jackHWE
#' @method plot jack.influential
#' @export
#' 
plot.jack.influential <- function(x, main = "", ...) {
  or <- as.vector(x$odds.ratio)
  or <- or[!is.na(or) | !is.nan(or) | !is.infinite(or)]
  odds.sort <- sort(or)
  odds.freq <- 1:length(odds.sort) / length(odds.sort)
  
  plot(odds.sort, odds.freq, type = "l", bty = "l", main = main,
       xlab = "Odds-ratio", ylab = "Cumulative frequency", ...
  )
  if(!all(is.null(x$influential))) {
    if(nrow(x$influential) > 0) {
      min.influential <- min(x$influential$odds.ratio)
      min.freq <- length(or[or < min.influential]) / length(or)
      abline(v = min.influential, lty = "dashed", lwd = 1)    
      min.txt <- paste(round(min.influential, 2), " = ", 
                       round(100 * min.freq, 0), "%", sep = "")
      text(min.influential, min.freq, min.txt, adj = c(1.1, 1.4), cex = 1)        
    }
  }
}